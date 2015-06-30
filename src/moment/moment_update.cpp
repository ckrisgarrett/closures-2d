/*
 File:   update.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "moment_solver.h"
#include <omp.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>
#include "../utils.h"


/*
    Scales moments on the entire grid.
*/
void MomentSolver::scaleMoments(double *moments)
{
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            for(int k = 0; k < c_numMoments; k++)
            {
                moments[I3D(i,j,k)] = moments[I3D(i,j,k)] * c_momentFilter[k];
            }
        }
    }
}


/*
    Scales flux on the entire grid.
*/
void MomentSolver::solveFlux(double *moments, double *flux, double dx, double dy)
{
    static bool firstTime = true;
    static double *tempVectors;
    static double *eastWestFluxVectors;
    static double *northSouthFluxVectors;
    if(firstTime)
    {
        firstTime = false;
        tempVectors           = (double*)malloc(c_numMoments * c_numOmpThreads * sizeof(double));
        eastWestFluxVectors   = (double*)malloc(c_numMoments * c_numOmpThreads * sizeof(double));
        northSouthFluxVectors = (double*)malloc(c_numMoments * c_numOmpThreads * sizeof(double));
    }
    
    
    int numX = c_gX[3] - c_gX[0] + 1;
    int numY = c_gY[3] - c_gY[0] + 1;
    int numGridPoints = numX * numY;
    
    
    // Go through the grid in one index (makes openmp faster).
    #pragma omp parallel for schedule(dynamic,1)
    for(int index = 0; index < numGridPoints; index++)
    {
        // Get 2d index.
        int i = index / numY;
        int j = index % numY;
        if(i < c_gX[1] || i > c_gX[2] || j < c_gY[1] || j > c_gY[2])
            continue;

        
        // Get temporary data for this thread.
        int threadNum = 0;
        #ifdef USE_OPENMP
        threadNum = omp_get_thread_num();
        #endif
        double *temp           = &tempVectors[threadNum * c_numMoments];
        double *eastWestFlux   = &eastWestFluxVectors[threadNum * c_numMoments];
        double *northSouthFlux = &northSouthFluxVectors[threadNum * c_numMoments];
        char tChar = 'T';
        double one = 1.0;
        double zero = 0.0;
        int incX = 1;
        
        // East West Flux
        for(int k = 0; k < c_numMoments; k++)
        {
            temp[k] = 0.25*moments[I3D(i+1,j,k)] + 0.75*moments[I3D(i,j,k)] - 
                1.25*moments[I3D(i-1,j,k)] + 0.25*moments[I3D(i-2,j,k)];
            //temp[k] = -0.125*moments[I3D(i+2,j,k)] + 0.75*moments[I3D(i+1,j,k)] - 
              //  0.75*moments[I3D(i-1,j,k)] + 0.125*moments[I3D(i-2,j,k)];
        }
        dgemv_(&tChar, &c_numMoments, &c_numMoments, &one, 
            c_pnFluxOperatorXiPlus, &c_numMoments, temp, &incX, &zero, eastWestFlux, &incX);
        
        for(int k = 0; k < c_numMoments; k++)
        {
            temp[k] = -0.25*moments[I3D(i+2,j,k)] + 1.25*moments[I3D(i+1,j,k)] - 
                0.75*moments[I3D(i,j,k)] - 0.25*moments[I3D(i-1,j,k)];
            //temp[k] = -0.125*moments[I3D(i+2,j,k)] + 0.75*moments[I3D(i+1,j,k)] - 
              //  0.75*moments[I3D(i-1,j,k)] + 0.125*moments[I3D(i-2,j,k)];
        }
        dgemv_(&tChar, &c_numMoments, &c_numMoments, &one, 
            c_pnFluxOperatorXiMinus, &c_numMoments, temp, &incX, &one, eastWestFlux, &incX);

        
        // North South Flux
        for(int k = 0; k < c_numMoments; k++)
        {
            temp[k] = 0.25*moments[I3D(i,j+1,k)] + 0.75*moments[I3D(i,j,k)] - 
                1.25*moments[I3D(i,j-1,k)] + 0.25*moments[I3D(i,j-2,k)];
            //temp[k] = -0.125*moments[I3D(i,j+2,k)] + 0.75*moments[I3D(i,j+1,k)] - 
              //  0.75*moments[I3D(i,j-1,k)] + 0.125*moments[I3D(i,j-2,k)];
        }
        dgemv_(&tChar, &c_numMoments, &c_numMoments, &one, 
            c_pnFluxOperatorEtaPlus, &c_numMoments, temp, &incX, &zero, northSouthFlux, &incX);
        
        for(int k = 0; k < c_numMoments; k++)
        {
            temp[k] = -0.25*moments[I3D(i,j+2,k)] + 1.25*moments[I3D(i,j+1,k)] - 
                0.75*moments[I3D(i,j,k)] - 0.25*moments[I3D(i,j-1,k)];
            //temp[k] = -0.125*moments[I3D(i,j+2,k)] + 0.75*moments[I3D(i,j+1,k)] - 
              //  0.75*moments[I3D(i,j-1,k)] + 0.125*moments[I3D(i,j-2,k)];
        }
        dgemv_(&tChar, &c_numMoments, &c_numMoments, &one, 
            c_pnFluxOperatorEtaMinus, &c_numMoments, temp, &incX, &one, northSouthFlux, &incX);
        

        // Add fluxes together.
        for(int k = 0; k < c_numMoments; k++)
        {
            flux[I3D(i,j,k)] = 1.0 / dx * eastWestFlux[k] + 1.0 / dy * northSouthFlux[k];
        }
    }
}


/*
    Updates one time step.
*/
void MomentSolver::update(double dt, double dx, double dy)
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
    if(c_filterType != FILTER_TYPE_NONE)
        scaleMoments(c_moments);
    communicateBoundaries();
    solveFlux(c_moments, c_flux, dx, dy);

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


    // Do first Euler step.
    if(c_filterType != FILTER_TYPE_NONE)
        scaleMoments(c_moments);
    communicateBoundaries();
    solveFlux(c_moments, c_flux, dx, dy);

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



void MomentSolver::periodicBoundary(double *u)
{
    int numX = (c_gX[3]-c_gX[0]+1);
    int numY = (c_gY[3]-c_gY[0]+1);
    int numGC = c_gX[1] - c_gX[0];
    
    for(int gc = 0; gc < numGC; gc++) {
    for(int j = 0; j < numY; j++) {
    for(int k = 0; k < c_numMoments; k++) {
        u[I3D(gc, j, k)] = u[I3D(c_gX[2] - numGC + 1 + gc, j, k)];
        u[I3D(c_gX[2] + 1 + gc, j, k)] = u[I3D(c_gX[1] + gc, j, k)];
    }}}
    
    for(int gc = 0; gc < numGC; gc++) {
    for(int i = 0; i < numX; i++) {
    for(int k = 0; k < c_numMoments; k++) {
        u[I3D(i, gc, k)] = u[I3D(i, c_gY[2] - numGC + 1 + gc, k)];
        u[I3D(i, c_gY[2] + 1 + gc, k)] = u[I3D(i, c_gY[1] + gc, k)];
    }}}
}


/*
    Updates one time step.
*//*
void MomentSolver::update(double dt, double dx, double dy)
{
    static bool firstTime = true;
    static double *u0;
    static double *u1;
    static double *u2;
    static double *u3;
    static double *flux1;
    static double *flux2;
    static double *flux3;
    if(firstTime)
    {
        firstTime = false;
        u0 = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numMoments * sizeof(double));
        u1 = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numMoments * sizeof(double));
        u2 = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numMoments * sizeof(double));
        u3 = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numMoments * sizeof(double));
        flux1 = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numMoments * sizeof(double));
        flux2 = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numMoments * sizeof(double));
        flux3 = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numMoments * sizeof(double));
    }
    
    
    // u0
    for(int i = c_gX[1]; i <= c_gX[2]; i++) {
    for(int j = c_gY[1]; j <= c_gY[2]; j++) {
    for(int k = 0; k < c_numMoments; k++) {
        u0[I3D(i,j,k)] = c_moments[I3D(i,j,k)];
    }}}
    periodicBoundary(u0);
    
    
    // Do first step.
    communicateBoundaries();
    for(int i = c_gX[1]; i <= c_gX[2]; i++) {
    for(int j = c_gY[1]; j <= c_gY[2]; j++) {
        for(int k = 0; k < c_numMoments; k++) {
            u1[I3D(i,j,k)] = u0[I3D(i,j,k)];
        }
        
        u1[I3D(i,j,0)] = u1[I3D(i,j,0)] / (1.0 - dt / 4.0 * (c_sigmaS[I2D(i,j)] - c_sigmaT[I2D(i,j)]));
        for(int k = 1; k < c_numMoments; k++) {
            u1[I3D(i,j,k)] = u1[I3D(i,j,k)] / (1.0 - dt / 4.0 * (-c_sigmaT[I2D(i,j)]));
        }
    }}
    periodicBoundary(u1);
    
    
    // Do second step.
    communicateBoundaries();
    solveFlux(u1, flux1, dx, dy);
    for(int i = c_gX[1]; i <= c_gX[2]; i++) {
    for(int j = c_gY[1]; j <= c_gY[2]; j++) {
        for(int k = 0; k < c_numMoments; k++) {
            u2[I3D(i,j,k)] = u0[I3D(i,j,k)] - dt / 2.0 * flux1[I3D(i,j,k)];
        }
        
        u2[I3D(i,j,0)] = u2[I3D(i,j,0)] / (1.0 - dt / 4.0 * (c_sigmaS[I2D(i,j)] - c_sigmaT[I2D(i,j)]));
        for(int k = 1; k < c_numMoments; k++) {
            u2[I3D(i,j,k)] = u2[I3D(i,j,k)] / (1.0 - dt / 4.0 * (-c_sigmaT[I2D(i,j)]));
        }
    }}
    periodicBoundary(u2);


    // Do third step.
    communicateBoundaries();
    solveFlux(u2, flux2, dx, dy);
    for(int i = c_gX[1]; i <= c_gX[2]; i++) {
    for(int j = c_gY[1]; j <= c_gY[2]; j++) {
        for(int k = 0; k < c_numMoments; k++) {
            u3[I3D(i,j,k)] = u0[I3D(i,j,k)] 
                - dt / 2.0 * flux1[I3D(i,j,k)] 
                - dt / 2.0 * flux2[I3D(i,j,k)] 
                - dt / 3.0 * c_sigmaT[I2D(i,j)] * u1[I3D(i,j,k)] 
                - dt / 3.0 * c_sigmaT[I2D(i,j)] * u2[I3D(i,j,k)];
        }
        u3[I3D(i,j,0)] += dt / 3.0 * c_sigmaS[I2D(i,j)] * u1[I3D(i,j,0)] 
                        + dt / 3.0 * c_sigmaS[I2D(i,j)] * u2[I3D(i,j,0)];
        
        u3[I3D(i,j,0)] = u3[I3D(i,j,0)] / (1.0 - dt / 3.0 * (c_sigmaS[I2D(i,j)] - c_sigmaT[I2D(i,j)]));
        for(int k = 1; k < c_numMoments; k++) {
            u3[I3D(i,j,k)] = u3[I3D(i,j,k)] / (1.0 - dt / 3.0 * (-c_sigmaT[I2D(i,j)]));
        }
    }}
    periodicBoundary(u3);
    
    
    // Do final step.
    communicateBoundaries();
    solveFlux(u3, flux3, dx, dy);
    for(int i = c_gX[1]; i <= c_gX[2]; i++) {
    for(int j = c_gY[1]; j <= c_gY[2]; j++) {
        for(int k = 0; k < c_numMoments; k++) {
            c_moments[I3D(i,j,k)] = u0[I3D(i,j,k)] 
                - dt / 3.0 * flux1[I3D(i,j,k)] 
                - dt / 3.0 * flux2[I3D(i,j,k)] 
                - dt / 3.0 * flux3[I3D(i,j,k)] 
                - dt / 3.0 * c_sigmaT[I2D(i,j)] * u1[I3D(i,j,k)] 
                - dt / 3.0 * c_sigmaT[I2D(i,j)] * u2[I3D(i,j,k)] 
                - dt / 3.0 * c_sigmaT[I2D(i,j)] * u3[I3D(i,j,k)];
        }
        c_moments[I3D(i,j,0)] += dt / 3.0 * c_sigmaS[I2D(i,j)] * u1[I3D(i,j,0)] 
                               + dt / 3.0 * c_sigmaS[I2D(i,j)] * u2[I3D(i,j,0)]
                               + dt / 3.0 * c_sigmaS[I2D(i,j)] * u3[I3D(i,j,0)];
    }}
    periodicBoundary(c_moments);
    
    
    if(c_filterType != FILTER_TYPE_NONE)
        scaleMoments(c_moments);
}*/

