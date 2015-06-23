/*
 File:   fobj_cuda.cu
 Author: Kris Garrett
 Date:   June 7, 2012
 
 Mimics fobj.cpp
 
 h CANNOT BE NULL
 This is important.  fobj should call this routine only when
 H != NULL
 
 g can be NULL
*/

#include "fobj_cuda.h"
#include <stdio.h>
#include "../../utils.h"


// The device variables.
// These are allocated on the CUDA cards in fobj_initialize.
static double **wDev;
static double **pDev;
static double **alphaDev;
static double **hDev;
static double **tempDev;

// Pinned memory and streams to use CUDA Streams.
static double **alphaPinned;
static double **hPinned;
static double **pPinned;
static cudaStream_t *stream;

// Used to determine if CUDA should be used.
static int numGpu;
static int numThreadsPerGpu;


static void checkError(cudaError_t error, int line)
{
    if(error != cudaSuccess)
        printf("%d: %s\n", line, cudaGetErrorString(error));
}


/*
    Fills in the array e^(alpha^T p).
    The array is of size quadrature points.
*/
__global__
void kernel_fillTempArrayMn(int nm, int nq, double *alpha, double *w, double *p,
                            double *tempArray)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if(k < nq)
    {
        tempArray[k] = 0;
        for(int i = 0; i < nm; i++)
        {
            tempArray[k] = tempArray[k] + alpha[i] * p[k * nm + i];
        }
        tempArray[k] = exp(tempArray[k]) * w[k];
    }
}


/*
    Calculates the Hessian for the Mn algorithm.
*/
__global__
void kernel_setHMn(int nm, int nq, double *tempArray, double *p, 
                   double *h)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int i = index / nm;
    int j = index % nm;
    if(index < nm * nm)
    {
        h[i * nm + j] = 0;
        
        // Do the integral part of the calculation.
        for(int k = 0; k < nq; k++)
        {
            h[i * nm + j] = h[i * nm + j] + 
                tempArray[k] * p[k * nm + i] * p[k * nm + j];
        }
    }
}


/*
    Calculates alpha^T p.
*/
__global__
void kernel_fillTempArrayPPn(int nm, int nq, double *alpha, double *p,
                             double *tempArray)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if(k < nq)
    {
        tempArray[k] = 0;
        for(int i = 0; i < nm; i++)
        {
            tempArray[k] = tempArray[k] + alpha[i] * p[k * nm + i];
        }
    }
}


/*
    Calculates the Hessian for the PPn algorithm.
*/
__global__
void kernel_setHPPn(int nm, int nq, double *tempArray, double *w, double *p, double delta, 
                    double *h)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int i = index / nm;
    int j = index % nm;
    if(index < nm * nm)
    {
        h[i * nm + j] = 0;
        
        // Do the integral part of the calculation.
        for(int k = 0; k < nq; k++)
        {
            double t_k = tempArray[k];
            h[i * nm + j] = h[i * nm + j] + 0.5 * w[k] * p[k * nm + i] * p[k * nm + j] * 
                           (1.0 + t_k / sqrt(t_k * t_k + 4.0 * delta));
        }
    }
}


/*
    Determines if CUDA should be used.
*/
extern "C"
int fobjUseCuda(int thread)
{
    // Return true
    if(thread < numGpu * numThreadsPerGpu)
        return 1;

    // Return false
    return 0;
}


/*
    Initialization that must be done in serial.
*/
extern "C"
void fobjCudaInitSerial(int numCudaCards, int numThreadsPerCudaCard)
{
    // Check for enough cuda cards and allocate memory for arrays.
    int numCards;
    checkError(cudaGetDeviceCount(&numCards), __LINE__);
    if(numCards < numCudaCards)
    {
        printf("fobj_cuda: Not enough CUDA cards in machine.\n");
        printf("fobj_cuda: Number of Cards: %d\n", numCards);
        printf("fobj_cuda: numCudaCards defined in the input deck: %d\n", numCudaCards);
        utils_abort();
    }

    wDev = (double**)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(double*));
    pDev = (double**)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(double*));
    alphaDev = (double**)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(double*));
    hDev = (double**)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(double*));
    tempDev = (double**)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(double*));
    alphaPinned = (double**)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(double*));
    hPinned = (double**)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(double*));
    pPinned = (double**)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(double*));
    stream = (cudaStream_t*)malloc(numCudaCards * numThreadsPerCudaCard * sizeof(cudaStream_t));
    
    numGpu = numCudaCards;
    numThreadsPerGpu = numThreadsPerCudaCard;
}




/*
    Initialization that must be done in parallel for openmp.
*/
extern "C"
void fobjCudaInitParallel(int nm, int nq, double *w, int thread)
{
    // Allocate memory on card and set quadrature.
    if(thread < numGpu * numThreadsPerGpu)
    {
        int gpu = thread % numGpu;
        checkError(cudaSetDevice(gpu), __LINE__);
        checkError(cudaStreamCreate(&stream[thread]), __LINE__);

        checkError(cudaMalloc((void**)&wDev[thread], nq * sizeof(double)), __LINE__);
        checkError(cudaMalloc((void**)&pDev[thread], nq * nm * sizeof(double)), __LINE__);
        checkError(cudaMalloc((void**)&alphaDev[thread], nm * sizeof(double)), __LINE__);
        checkError(cudaMalloc((void**)&hDev[thread],     nm * nm * sizeof(double)), __LINE__);
        checkError(cudaMalloc((void**)&tempDev[thread],  nq * sizeof(double)), __LINE__);
        
        checkError(cudaMemcpy(wDev[thread], w, nq * sizeof(double), cudaMemcpyHostToDevice), __LINE__);

        checkError(cudaHostAlloc((void**)&alphaPinned[thread], nm * sizeof(double), 
            cudaHostAllocPortable), __LINE__);
        checkError(cudaHostAlloc((void**)&hPinned[thread], nm * nm * sizeof(double), 
            cudaHostAllocPortable), __LINE__);
        checkError(cudaHostAlloc((void**)&pPinned[thread], nq * nm * sizeof(double), 
            cudaHostAllocPortable), __LINE__);
    }
}


/*
    fobj replacement using CUDA for the Mn algorithm.
*/
extern "C"
double fobjMnCuda(int nm, int nq, double *alpha, double *u, double *p, int thread, 
                  double *g, double *h)
{
    int numThreadsPerBlock = 256;
    int numBlocks = (nq-1) / numThreadsPerBlock + 1;
    int numThreadsPerBlockH = 256;
    int numBlocksH = (nm*nm-1) / numThreadsPerBlockH + 1;
    
    // Set alpha on card and solve for H.
    memcpy(alphaPinned[thread], alpha, nm * sizeof(double));
    memcpy(pPinned[thread], p, nm * nq * sizeof(double));
    checkError(cudaMemcpyAsync(alphaDev[thread], alphaPinned[thread], nm * sizeof(double), 
               cudaMemcpyHostToDevice, stream[thread]), __LINE__);
    checkError(cudaMemcpyAsync(pDev[thread], pPinned[thread], nm * nq * sizeof(double), 
               cudaMemcpyHostToDevice, stream[thread]), __LINE__);
    
    kernel_fillTempArrayMn<<<numBlocks, numThreadsPerBlock, 0, stream[thread]>>>(nm, nq, 
                          alphaDev[thread], wDev[thread], pDev[thread], 
                          tempDev[thread]);
    kernel_setHMn<<<numBlocksH, numThreadsPerBlockH, 0, stream[thread]>>>(nm, nq, 
                 tempDev[thread], pDev[thread], hDev[thread]);
    
    checkError(cudaMemcpyAsync(hPinned[thread], hDev[thread], nm * nm * sizeof(double), 
               cudaMemcpyDeviceToHost, stream[thread]), __LINE__);
    checkError(cudaStreamSynchronize(stream[thread]), __LINE__);
    memcpy(h, hPinned[thread], nm * nm * sizeof(double));

    
    // Get g.
    double p00 = p[0];
    if(g != NULL);
    {
        for(int k = 0; k < nm; k++)
            g[k] = h[k] / p00 - u[k];
    }
    
    // Get f.
    double f = h[0] / (p00 * p00);
    
    // f = f - alpha^T u
    for(int k = 0; k < nm; k++)
        f = f - alpha[k] * u[k];

    return f;
}


/*
    Calculates f and g for the PPn algorithm.
*/
static
double getfgPPn(int nm, int nq, double *alpha, double *u, double *w, double *p, 
                double delta, double *g)
{
    double f = 0;
    for(int k = 0; k < nm; k++)
        g[k] = 0.0;

    // Integrate.
    for(int q = 0; q < nq; q++)
    {
        // Compute temp = alpha^T * p.
        double temp = 0;
        for(int i = 0; i < nm; i++)
        {
            temp = temp + alpha[i] * p[q*nm+i];
        }

        // Compute temp2 = sqrt((alpha^t * p)^2 + 4 * delta).
        double temp2 = sqrt(temp * temp + 4 * delta);

        // Compute temp3 = 0.5 alpha^T * p + 0.5 sqrt((alpha^T * p)^2 + 4 * delta).
        double temp3 = 0.5 * (temp + temp2);

        // Compute f.
        f = f + w[q] * (0.5 * temp3 * temp3 - delta + delta * log(temp3));

        // Compute g.
        double temp4 = w[q] * temp3;
        for(int i = 0; i < nm; i++)
            g[i] = g[i] + p[q*nq+i] * temp4;
    }

    // Finish the computations not requiring integration.
    for(int i = 0; i < nm; i++)
    {
        f = f - alpha[i] * u[i];

        if(g != NULL)
            g[i] = g[i] - u[i];
    }

    return f;
}


/*
    fobj replacement using CUDA for the PPn algorithm.
*/
extern "C"
double fobjPPnCuda(int nm, int nq, double *alpha, double *u, double *w, double *p, 
                   double delta, int thread, double *g, double *h)
{
    int numThreadsPerBlock = 96;
    int numBlocks = (nq-1) / numThreadsPerBlock + 1;
    int numThreadsPerBlockH = 96;
    int numBlocksH = (nm*nm-1) / numThreadsPerBlockH + 1;
    
    // Set alpha on card and solve for H.
    memcpy(alphaPinned[thread], alpha, nm * sizeof(double));
    memcpy(pPinned[thread], p, nm * nq * sizeof(double));
    cudaMemcpyAsync(alphaDev[thread], alphaPinned[thread], nm * sizeof(double), 
               cudaMemcpyHostToDevice, stream[thread]);
    cudaMemcpyAsync(pDev[thread], pPinned[thread], nm * nq * sizeof(double), 
               cudaMemcpyHostToDevice, stream[thread]);
    
    kernel_fillTempArrayPPn<<<numBlocks, numThreadsPerBlock, 0, stream[thread]>>>
        (nm, nq, alphaDev[thread], pDev[thread], tempDev[thread]);
    kernel_setHPPn<<<numBlocksH, numThreadsPerBlockH, 0, stream[thread]>>>
        (nm, nq, tempDev[thread], wDev[thread], pDev[thread], delta, hDev[thread]);
    
    double f = getfgPPn(nm, nq, alpha, u, w, p, delta, g);
    
    cudaMemcpyAsync(hPinned[thread], hDev[thread], nm * nm * sizeof(double), 
               cudaMemcpyDeviceToHost, stream[thread]);
    cudaStreamSynchronize(stream[thread]);
    memcpy(h, hPinned[thread], nm * nm * sizeof(double));
    

    return f;
}

