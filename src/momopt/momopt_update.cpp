/*
 File:   update.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "momopt_solver.h"
#include "../utils.h"
#include <gsl/gsl_blas.h>
#include <omp.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef USE_CUDA_FLUX
#include "solve_flux_cuda.h"
#endif


// Temporary variables.
static double *ansatzGrid;


/*
    Initialize temporary variables.
*/
void MomOptSolver::initUpdate(int numMoments)
{
    ansatzGrid = new double[(c_gX[3] - c_gX[0] + 1) * (c_gY[3] - c_gY[0] + 1)];
}


/* 
    Given alpha and the quadrature point, calculate the ansatz at that quadrature point.
*/
double MomOptSolver::computeAnsatz(double *alpha, int q)
{
    double kinetic = 0.0;
    for(int k = 0; k < c_numMoments; k++)
    {
        kinetic += alpha[k] * c_spHarm[q*c_numManyMoments+k];
    }
    
    switch(c_momentType)
    {
        case MOMENT_TYPE_MN:
            return exp(kinetic);
        case MOMENT_TYPE_PPN:
            if(kinetic > 0)
                return 0.5 * kinetic + 0.5 * sqrt(kinetic * kinetic + 4 * c_deltaPPn);
            return -c_deltaPPn / (0.5 * kinetic - 0.5 * sqrt(kinetic*kinetic + 4*c_deltaPPn));
        default:
            printf("computeAnsatz: momentType not in range\n");
            utils_abort();
    }

    // Should never reach this point.
    return 0.0;
}


/*
    Returns 0 if xy < 1
    Returns min(x,y) if x > 0 and y > 0
    Returns -min(|x|,|y|) if x < 0 and y < 0
*/
static double minmod(double x, double y)
{
    return SIGN(1.0,x)* MAX(0.0, MIN(fabs(x),y*SIGN(1.0,x) ) );
}


/*
    Computes the undivided differences with limiter.
*/
static double slopefit(double left, double center, double right, double theta)
{
    return minmod(theta*(right-center),
                  minmod(0.5*(right-left), theta*(center-left)) );
}


/*
    Solves for the flux at each grid cell.
*/
void MomOptSolver::solveFlux(double *moments, double *flux, double *alpha, double dx, double dy)
{
    int numX = c_gX[3] - c_gX[0] + 1;
    int numY = c_gY[3] - c_gY[0] + 1;
    int numGridPoints = numX * numY;
    
    
    // If using CUDA ...
    #ifdef USE_CUDA_FLUX
    // Solve for the flux on the GPU.
    solveFlux_cuda(numX, numY, c_numMoments, c_numQuadPoints, c_theta, dx, dy, alpha, flux);
    
    
    // If not using CUDA ...
    #else
    // Zero flux.
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            for(int k = 0; k < c_numMoments; k++)
                flux[I3D(i,j,k)] = 0.0;
        }
    }

    // Add up each quadrature's flux.
    for(int q = 0; q < c_numQuadPoints; q++)
    {
        #pragma omp parallel for schedule(dynamic,1)
        for(int index = 0; index < numGridPoints; index++)
        {
            int i = index / numY;
            int j = index % numY;
            
            ansatzGrid[index] = computeAnsatz(&alpha[I3D(i,j,0)], q);
        }
        
        #pragma omp parallel for schedule(dynamic,1)
        for(int index = 0; index < numGridPoints; index++)
        {
            int i = index / numY;
            int j = index % numY;
            if(i < c_gX[1] || i > c_gX[2] || j < c_gY[1] || j > c_gY[2])
                continue;
            
            
            double xi  = c_xi[q];
            double eta = c_eta[q];
            double w   = c_quadWeights[q];
            double flux1 = 0;
            double flux2 = 0;
            double flux3 = 0;
            double flux4 = 0;
            
            if(xi > 0)
            {
                double k1 = ansatzGrid[(i-2) * numY + j];
                double k2 = ansatzGrid[(i-1) * numY + j];
                double k3 = ansatzGrid[(i)   * numY + j];
                double k4 = ansatzGrid[(i+1) * numY + j];
                
                flux1 = k3 + 0.5 * slopefit(k2, k3, k4, c_theta) - 
                    k2 - 0.5 * slopefit(k1, k2, k3, c_theta);
                flux1 = flux1 * xi * w;
            }
            if(xi < 0)
            {
                double k1 = ansatzGrid[(i-1) * numY + j];
                double k2 = ansatzGrid[(i)   * numY + j];
                double k3 = ansatzGrid[(i+1) * numY + j];
                double k4 = ansatzGrid[(i+2) * numY + j];
                
                flux2 = k3 - 0.5 * slopefit(k2, k3, k4, c_theta) - 
                    k2 + 0.5 * slopefit(k1, k2, k3, c_theta);
                flux2 = flux2 * xi * w;
            }
            if(eta > 0)
            {
                double k1 = ansatzGrid[i * numY + (j-2)];
                double k2 = ansatzGrid[i * numY + (j-1)];
                double k3 = ansatzGrid[i * numY + (j)];
                double k4 = ansatzGrid[i * numY + (j+1)];
                
                flux3 = k3 + 0.5 * slopefit(k2, k3, k4, c_theta) - 
                    k2 - 0.5 * slopefit(k1, k2, k3, c_theta);
                flux3 = flux3 * eta * w;
            }
            if(eta < 0)
            {
                double k1 = ansatzGrid[i * numY + (j-1)];
                double k2 = ansatzGrid[i * numY + (j)];
                double k3 = ansatzGrid[i * numY + (j+1)];
                double k4 = ansatzGrid[i * numY + (j+2)];
                
                flux4 = k3 - 0.5 * slopefit(k2, k3, k4, c_theta) - 
                    k2 + 0.5 * slopefit(k1, k2, k3, c_theta);
                flux4 = flux4 * eta * w;
            }
            
            
            for(int k = 0; k < c_numMoments; k++)
            {
                flux[I3D(i,j,k)] += c_spHarm[q*c_numManyMoments+k] * 
                    ((flux1 + flux2) / dx + (flux3 + flux4) / dy);
            }
        }
    }
    #endif
}


/*
    Solves the optimization problem on the entire grid.
*/
void MomOptSolver::solveOptimization()
{
    int numX = c_gX[3] - c_gX[0] + 1;
    int numY = c_gY[3] - c_gY[0] + 1;
    int numGridPoints = numX * numY;
    
    #pragma omp parallel for schedule(dynamic,1)
    for(int index = 0; index < numGridPoints; index++)
    {
        int i = index / numY;
        int j = index % numY;

        OPTIONS options;
        options.momentType = c_momentType;
        options.maxIter = c_maxIter;
        options.maxBfgsIter = c_maxIterBfgs;
        options.tolAbs = c_tolerance * c_tolerance;
        options.tolRel = c_tolerance;
        options.tolGamma = 1.1;
        options.condHMax = c_condHMax;
        options.condHMaxBfgs = c_condHMaxBfgs;
        options.deltaPPn = c_deltaPPn;
        options.useClebschGordan = (c_useClebschGordan == 1) ? true : false;
        
        OUTS outs;
        switch(c_optType)
        {
            case OPT_TYPE_ORIGINAL:
                opt(c_numMoments, c_numManyMoments, c_numQuadPoints, 
                    &c_moments[I3D(i,j,0)], &c_alpha[I3D(i,j,0)], 
                    options, c_quadWeights, c_spHarm, &outs);
                break;
            case OPT_TYPE_CHANGE_BASIS:
                optcb(c_numMoments, c_numQuadPoints, &c_moments[I3D(i,j,0)], 
                      &c_P2M[I3D2(i,j,0)], &c_alphaP[I3D(i,j,0)], &c_alpha[I3D(i,j,0)], 
                      c_quadWeights, c_spHarm, options, &outs);
                break;
            case OPT_TYPE_BFGS:
                optbfgs(c_numMoments, c_numQuadPoints, &c_moments[I3D(i,j,0)], 
                        &c_P2M[I3D2(i,j,0)], &c_alphaP[I3D(i,j,0)], &c_alpha[I3D(i,j,0)],  
                        c_quadWeights, c_spHarm, options, &outs);
                break;
            default:
                printf("update.cpp: optType out of bounds.\n");
                utils_abort();
        }
        
        
        #pragma omp critical
        {
            if(outs.iter > c_optStats.maxIter)
                c_optStats.maxIter = outs.iter;
            if(outs.iterGamma > c_optStats.maxGammaIter)
                c_optStats.maxGammaIter = outs.iterGamma;
            for(int r = 0; r < NUM_REGULARIZATIONS; r++)
            {
                if(fabs(outs.r - REGULARIZATIONS[r]) < 1e-10)
                    c_optStats.histReg[r]++;
            }
            
            c_optStats.iterMean = 
                (c_optStats.numDualSolves * c_optStats.iterMean + outs.iter) / 
                (c_optStats.numDualSolves + 1.0);
            c_optStats.iterGammaMean = 
                (c_optStats.numDualSolves * c_optStats.iterGammaMean + outs.iterGamma) / 
                (c_optStats.numDualSolves + 1.0);
            c_optStats.numDualSolves++;
        }
    }
}


/*
    Updates one time step.
*/
void MomOptSolver::update(double dt, double dx, double dy)
{
    static bool firstTime = true;
    static double *momentsOld;
    if(firstTime)
    {
        firstTime = false;
        momentsOld = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numMoments * sizeof(double));
    }
    
    
    // Copy initial grid for use later.
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            for(int k = 0; k < c_numMoments; k++)
            {
                momentsOld[I3D(i,j,k)] = c_moments[I3D(i,j,k)];
            }
        }
    }
    

    // Do first Euler step.
    communicateBoundaries();
    solveOptimization();
    solveFlux(c_moments, c_flux, c_alpha, dx, dy);

    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            c_moments[I3D(i,j,0)] = c_moments[I3D(i,j,0)] * 
                (1.0 + dt * c_sigmaS[I2D(i,j)] - dt * c_sigmaT[I2D(i,j)]) - dt * c_flux[I3D(i,j,0)];
            
            for(int k = 1; k < c_numMoments; k++)
            {
                c_moments[I3D(i,j,k)] = c_moments[I3D(i,j,k)] * 
                    (1.0 - dt * c_sigmaT[I2D(i,j)]) - dt * c_flux[I3D(i,j,k)];
                if(c_initCond == 1) {
                    double x_i = c_initX + (i-c_gX[1]) * dx;
                    double y_j = c_initY + (j-c_gY[1]) * dy;
                    if(fabs(x_i) < 0.5 && fabs(y_j) < 0.5) {
                        c_moments[I3D(i,j,0)] = c_moments[I3D(i,j,0)] + dt * 2.0 * sqrt(M_PI);
                    }
                }
            }
        }
    }
    
    
    // Do second Euler step.
    communicateBoundaries();
    solveOptimization();
    solveFlux(c_moments, c_flux, c_alpha, dx, dy);
    
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            c_moments[I3D(i,j,0)] = c_moments[I3D(i,j,0)] * 
                (1.0 + dt * c_sigmaS[I2D(i,j)] - dt * c_sigmaT[I2D(i,j)]) - dt * c_flux[I3D(i,j,0)];
            
            for(int k = 1; k < c_numMoments; k++)
            {
                c_moments[I3D(i,j,k)] = c_moments[I3D(i,j,k)] * 
                    (1.0 - dt * c_sigmaT[I2D(i,j)]) - dt * c_flux[I3D(i,j,k)];
                if(c_initCond == 1) {
                    double x_i = c_initX + (i-c_gX[1]) * dx;
                    double y_j = c_initY + (j-c_gY[1]) * dy;
                    if(fabs(x_i) < 0.5 && fabs(y_j) < 0.5) {
                        c_moments[I3D(i,j,0)] = c_moments[I3D(i,j,0)] + dt * 2.0 * sqrt(M_PI);
                    }
                }
            }
        }
    }
    
    
    // Average initial with second Euler step.
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            for(int k = 0; k < c_numMoments; k++)
            {
                c_moments[I3D(i,j,k)] = 0.5 * (c_moments[I3D(i,j,k)] + momentsOld[I3D(i,j,k)]);
            }
        }
    }
}

