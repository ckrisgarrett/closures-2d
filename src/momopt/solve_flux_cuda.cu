/*
 File:   solve_flux_cuda.cu
 Author: Kris Garrett
 Date:   April 19, 2013
*/

#include "solve_flux_cuda.h"
#include "../utils.h"


// GPU device variables
static double *alphaDev;
static double *pDev;
static double *ansatzDev;
static double *wDev;
static double *xiDev;
static double *etaDev;
static double *fluxDev;
static double *fluxTempDev;


/*
    Initialization for solveFlux_cuda.

    nx:  size of grid in x-direction
    ny:  size of grid in y-direction
    nm:  number of moments
    w:   quadrature weights
    p:   moments at quadrature points
    xi:  xi values at quadrature points
    eta: eta values at quadrature points
*/
extern "C"
void solveFluxInit_cuda(int nx, int ny, int nm, int nm2, int nq, double *w, double *p, 
                        double *xi, double *eta)
{
    cudaMalloc((void**)&alphaDev,  nx * ny * nm * sizeof(double));
    cudaMalloc((void**)&pDev,      nq * nm * sizeof(double));
    cudaMalloc((void**)&ansatzDev, nx * ny * sizeof(double));
    cudaMalloc((void**)&wDev,      nq * sizeof(double));
    cudaMalloc((void**)&xiDev,     nq * sizeof(double));
    cudaMalloc((void**)&etaDev,    nq * sizeof(double));
    cudaMalloc((void**)&fluxDev,   nx * ny * nm * sizeof(double));
    cudaMalloc((void**)&fluxTempDev,   nx * ny * sizeof(double));
    
    // Must do copy one row at a time since nm2 may be bigger than number of moments
    // for Clebsch-Gordan.
    double *pSmall = (double*)malloc(nq * nm * sizeof(double));
    for(int q = 0; q < nq; q++)
    {
        memcpy(&pSmall[q * nm], &p[q * nm2], nm * sizeof(double));
    }
    cudaMemcpy(pDev,   pSmall,  nq * nm * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(wDev,   w,       nq * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(xiDev,  xi,      nq * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(etaDev, eta,     nq * sizeof(double), cudaMemcpyHostToDevice);
    
    free(pSmall);
}


/*
    Helper functions to calculate slopes.
*/
__device__
double minmod(double x, double y)
{
    return SIGN(1.0,x) * MAX(0.0, MIN(fabs(x), y * SIGN(1.0, x)));
}

__device__
double slopefit(double left, double center, double right, double theta)
{
    return minmod(theta*(right-center), minmod(0.5*(right-left), theta*(center-left)));
}


/*
    Computes the ansatz at quadrature index q for the entire grid.
    
    n:      total size of the grid
    nm:     number of moments
    q:      quadrature index
    alpha:  grid of alpha vectors
    p:      moments
    ansatz: grid of ansatzes
*/
/*__global__
void kernel_computeAnsatz(int n, int nm, int q, double *alpha, double *p, double *ansatz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(index < n)
    {
        double *alpha_i = &alpha[index * nm];
        double kinetic = 0.0;
        for(int k = 0; k < nm; k++)
        {
            kinetic += alpha_i[k] * p[q * nm + k];
        }
        ansatz[index] = exp(kinetic);
    }
}*/


__global__
void kernel_computeAnsatz(int n, int nm, int q, double *alpha, double *p, double *ansatz)
{
    __shared__ double temp[16][16];
    int space_index = blockIdx.x * blockDim.x + threadIdx.y;
    
    
    temp[threadIdx.y][threadIdx.x] = 0.0;
    __syncthreads();
    
    
    if(space_index < n)
    {
        for(int batch_index = 0; batch_index < nm; batch_index += 16)
        {
            int k = batch_index + threadIdx.x;
            
            if(k < nm)
                temp[threadIdx.y][threadIdx.x] += alpha[space_index * nm + k] * p[q*nm+k];
        }
        __syncthreads();
        
        if(threadIdx.y == 0)
        {
            double ans = 0.0;
            for(int i = 0; i < 16; i++)
            {
                ans += temp[threadIdx.x][i];
            }
            ans = exp(ans);
            ansatz[blockIdx.x * blockDim.x + threadIdx.x] = ans;
        }
    }
}


/*
    Solve flux on GPU.

    nx:     size of grid in x-direction
    ny:     size of grid in y-direction
    q:      quadrature index
    nm:     number of moments
    theta:  value for slopefit
    xi:     xi at quadrature points
    eta:    eta at quadrature points
    w:      quadrature weights
    p:      moments at quadrature points
    ansatz: grid of ansatzes
    flux:   grid of fluxes
*/
/*__global__
void kernel_solveFlux(int nx, int ny, int q, int nm, double theta, double dx, double dy, 
                      double *xi, double *eta, double *w, double *p, double *ansatz, double *flux)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int i = index / ny;
    int j = index % ny;
    
    if( (index < nx * ny) && (i > 2) && (i < nx - 2) && (j > 2) && (j < ny - 2) )
    {
        double flux1 = 0;
        double flux2 = 0;
        double xi_q = xi[q];
        double eta_q = eta[q];
        double w_q = w[q];
        
        if(xi_q > 0)
        {
            double k1 = ansatz[(i-2) * ny + j];
            double k2 = ansatz[(i-1) * ny + j];
            double k3 = ansatz[(i)   * ny + j];
            double k4 = ansatz[(i+1) * ny + j];
            
            flux1 = k3 + 0.5 * slopefit(k2, k3, k4, theta) - 
                k2 - 0.5 * slopefit(k1, k2, k3, theta);
            flux1 = flux1 * xi_q * w_q;
        }
        else
        {
            double k1 = ansatz[(i-1) * ny + j];
            double k2 = ansatz[(i)   * ny + j];
            double k3 = ansatz[(i+1) * ny + j];
            double k4 = ansatz[(i+2) * ny + j];
            
            flux1 = k3 - 0.5 * slopefit(k2, k3, k4, theta) - 
                k2 + 0.5 * slopefit(k1, k2, k3, theta);
            flux1 = flux1 * xi_q * w_q;
        }
        if(eta_q > 0)
        {
            double k1 = ansatz[i * ny + (j-2)];
            double k2 = ansatz[i * ny + (j-1)];
            double k3 = ansatz[i * ny + (j)];
            double k4 = ansatz[i * ny + (j+1)];
            
            flux2 = k3 + 0.5 * slopefit(k2, k3, k4, theta) - 
                k2 - 0.5 * slopefit(k1, k2, k3, theta);
            flux2 = flux2 * eta_q * w_q;
        }
        else
        {
            double k1 = ansatz[i * ny + (j-1)];
            double k2 = ansatz[i * ny + (j)];
            double k3 = ansatz[i * ny + (j+1)];
            double k4 = ansatz[i * ny + (j+2)];
            
            flux2 = k3 - 0.5 * slopefit(k2, k3, k4, theta) - 
                k2 + 0.5 * slopefit(k1, k2, k3, theta);
            flux2 = flux2 * eta_q * w_q;
        }
        
        
        for(int k = 0; k < nm; k++)
        {
            flux[index * nm + k] += p[q * nm + k] * (flux1 / dx + flux2 / dy);
        }
    }
}*/


__global__
void kernel_solveFluxTemp(int nx, int ny, int q, int nm, double theta, double dx, double dy, 
                      double *xi, double *eta, double *w, double *p, double *ansatz, double *fluxTemp)
{
    __shared__ double ansatz_shared[20][20];
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    
    
    ansatz_shared[threadIdx.y][threadIdx.x] = ansatz[(i-2) * ny + (j-2)];
    if(threadIdx.y < 4)
        ansatz_shared[threadIdx.y+16][threadIdx.x] = ansatz[(i-2+16) * ny + (j-2)];
    if(threadIdx.x < 4)
        ansatz_shared[threadIdx.y][threadIdx.x+16] = ansatz[(i-2) * ny + (j-2+16)];
    if(threadIdx.x < 4 && threadIdx.y < 4)
        ansatz_shared[threadIdx.y+16][threadIdx.x+16] = ansatz[(i-2+16) * ny + (j-2+16)];
    __syncthreads();
    
    
    if( (i > 2) && (i < nx - 2) && (j > 2) && (j < ny - 2) )
    {
        double flux1 = 0;
        double flux2 = 0;
        double xi_q = xi[q];
        double eta_q = eta[q];
        double w_q = w[q];
        
        if(xi_q > 0)
        {
            double k1 = ansatz_shared[threadIdx.y+2-2][threadIdx.x+2];
            double k2 = ansatz_shared[threadIdx.y+2-1][threadIdx.x+2];
            double k3 = ansatz_shared[threadIdx.y+2][threadIdx.x+2];
            double k4 = ansatz_shared[threadIdx.y+2+1][threadIdx.x+2];
            
            flux1 = k3 + 0.5 * slopefit(k2, k3, k4, theta) - 
                k2 - 0.5 * slopefit(k1, k2, k3, theta);
            flux1 = flux1 * xi_q * w_q;
        }
        else
        {
            double k1 = ansatz_shared[threadIdx.y+2-1][threadIdx.x+2];
            double k2 = ansatz_shared[threadIdx.y+2][threadIdx.x+2];
            double k3 = ansatz_shared[threadIdx.y+2+1][threadIdx.x+2];
            double k4 = ansatz_shared[threadIdx.y+2+2][threadIdx.x+2];
            
            flux1 = k3 - 0.5 * slopefit(k2, k3, k4, theta) - 
                k2 + 0.5 * slopefit(k1, k2, k3, theta);
            flux1 = flux1 * xi_q * w_q;
        }
        if(eta_q > 0)
        {
            double k1 = ansatz_shared[threadIdx.y+2][threadIdx.x+2-2];
            double k2 = ansatz_shared[threadIdx.y+2][threadIdx.x+2-1];
            double k3 = ansatz_shared[threadIdx.y+2][threadIdx.x+2];
            double k4 = ansatz_shared[threadIdx.y+2][threadIdx.x+2+1];
            
            flux2 = k3 + 0.5 * slopefit(k2, k3, k4, theta) - 
                k2 - 0.5 * slopefit(k1, k2, k3, theta);
            flux2 = flux2 * eta_q * w_q;
        }
        else
        {
            double k1 = ansatz_shared[threadIdx.y+2][threadIdx.x+2-1];
            double k2 = ansatz_shared[threadIdx.y+2][threadIdx.x+2];
            double k3 = ansatz_shared[threadIdx.y+2][threadIdx.x+2+1];
            double k4 = ansatz_shared[threadIdx.y+2][threadIdx.x+2+2];
            
            flux2 = k3 - 0.5 * slopefit(k2, k3, k4, theta) - 
                k2 + 0.5 * slopefit(k1, k2, k3, theta);
            flux2 = flux2 * eta_q * w_q;
        }
        
        fluxTemp[i*ny+j] = flux1 / dx + flux2 / dy;
        /*for(int k = 0; k < nm; k++)
        {
            flux[(i*ny+j) * nm + k] += p[q * nm + k] * (flux1 / dx + flux2 / dy);
        }*/
    }
}


__global__
void kernel_solveFlux(int nx, int ny, int nm, int q, double *fluxTemp, double *p, double *flux)
{
    int space_index = blockIdx.x * blockDim.x + threadIdx.y;
    int i = space_index / ny;
    int j = space_index % ny;
    
    if( (i > 2) && (i < nx - 2) && (j > 2) && (j < ny - 2) )
    {
        for(int batch_index = 0; batch_index < nm; batch_index += 16)
        {
            int k = threadIdx.x + batch_index;
            if(k < nm)
                flux[space_index * nm + k] += p[q * nm + k] * fluxTemp[space_index];
        }
    }
}


/*
    Solve flux on GPU driver.

    nx:    size of grid in x-direction
    ny:    size of grid in y-direction
    nm:    number of moments
    q:     quadrature index
    theta: value for slopefit
    alpha: grid of alpha vectors
    flux:  grid of fluxes
*/
extern "C"
void solveFlux_cuda(int nx, int ny, int nm, int nq, double theta, double dx, double dy, 
                    double *alpha, double *flux)
{
    int n = nx * ny;
    
    dim3 nt1(16,16);
    int nb1 = (n-1) / 16 + 1;
    dim3 nt2(16,16);
    dim3 nb2((ny-1)/16+1, (nx-1)/16+1);
    
    cudaMemcpy(alphaDev, alpha, n * nm * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(fluxDev, 0, n * nm * sizeof(double));
    for(int q = 0; q < nq; q++)
    {
        kernel_computeAnsatz<<<nb1, nt1>>>(n, nm, q, alphaDev, pDev, ansatzDev);
        kernel_solveFluxTemp<<<nb2, nt2>>>(nx, ny, q, nm, theta, dx, dy, 
                                       xiDev, etaDev, wDev, pDev, ansatzDev, fluxTempDev);
        kernel_solveFlux<<<nb1, nt1>>>(nx, ny, nm, q, fluxTempDev, pDev, fluxDev);
    }
    cudaMemcpy(flux, fluxDev, n * nm * sizeof(double), cudaMemcpyDeviceToHost);
}

