/*
 File:   update.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "dn_solver.h"
#include <omp.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>
#include "../utils.h"

#ifdef USE_PAPI
#include "../profiling.h"
#endif

/*
    Solves for the flux at each grid cell.
*/
void DnSolver::solveFlux(double *moments, double *flux, double dx, double dy)
{
    #ifdef USE_PAPI
    profile_start_update("DnSolver::solveFlux");
    #endif

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
        }
        dgemv_(&tChar, &c_numMoments, &c_numMoments, &one, 
            c_fluxOperatorXiPlus, &c_numMoments, temp, &incX, &zero, eastWestFlux, &incX);
        
        for(int k = 0; k < c_numMoments; k++)
        {
            temp[k] = -0.25*moments[I3D(i+2,j,k)] + 1.25*moments[I3D(i+1,j,k)] - 
                0.75*moments[I3D(i,j,k)] - 0.25*moments[I3D(i-1,j,k)];
        }
        dgemv_(&tChar, &c_numMoments, &c_numMoments, &one, 
            c_fluxOperatorXiMinus, &c_numMoments, temp, &incX, &one, eastWestFlux, &incX);

        
        // North South Flux
        for(int k = 0; k < c_numMoments; k++)
        {
            temp[k] = 0.25*moments[I3D(i,j+1,k)] + 0.75*moments[I3D(i,j,k)] - 
                1.25*moments[I3D(i,j-1,k)] + 0.25*moments[I3D(i,j-2,k)];
        }
        dgemv_(&tChar, &c_numMoments, &c_numMoments, &one, 
            c_fluxOperatorEtaPlus, &c_numMoments, temp, &incX, &zero, northSouthFlux, &incX);
        
        for(int k = 0; k < c_numMoments; k++)
        {
            temp[k] = -0.25*moments[I3D(i,j+2,k)] + 1.25*moments[I3D(i,j+1,k)] - 
                0.75*moments[I3D(i,j,k)] - 0.25*moments[I3D(i,j-1,k)];
        }
        dgemv_(&tChar, &c_numMoments, &c_numMoments, &one, 
            c_fluxOperatorEtaMinus, &c_numMoments, temp, &incX, &one, northSouthFlux, &incX);
        

        // Add fluxes together.
        for(int k = 0; k < c_numMoments; k++)
        {
            flux[I3D(i,j,k)] = 1.0 / dx * eastWestFlux[k] + 1.0 / dy * northSouthFlux[k];
        }
    }
    
    #pragma omp parallel for schedule(dynamic,1)
    for(int index = 0; index < numGridPoints; index++)
    {
        int i = index / numY;
        int j = index % numY;
        if(i < c_gX[1] || i > c_gX[2] || j < c_gY[1] || j > c_gY[2])
            continue;
        
        
        for(int k = 0; k < c_numMoments; k++)
        {
            double sigma_ip = (c_sigmaT[I2D(i+1,j)] + c_sigmaT[I2D(i,j)]) / 2.0;
            double sigma_in = (c_sigmaT[I2D(i,j)] + c_sigmaT[I2D(i-1,j)]) / 2.0;
            double sigma_jp = (c_sigmaT[I2D(i,j+1)] + c_sigmaT[I2D(i,j)]) / 2.0;
            double sigma_jn = (c_sigmaT[I2D(i,j)] + c_sigmaT[I2D(i,j-1)]) / 2.0;
            
            double dx_uip = 1.0 / sigma_ip * 
                (c_moments[I3D(i+1,j,k)] - c_moments[I3D(i,j,k)]) / dx;
            double dx_uin = 1.0 / sigma_in * 
                (c_moments[I3D(i,j,k)] - c_moments[I3D(i-1,j,k)]) / dx;
            double dy_uip = 1.0 / sigma_ip * 
                (c_moments[I3D(i+1,j+1,k)] + c_moments[I3D(i,j+1,k)] - 
                 c_moments[I3D(i+1,j-1,k)] - c_moments[I3D(i,j-1,k)]) / (4.0 * dy);
            double dy_uin = 1.0 / sigma_in * 
                (c_moments[I3D(i,j+1,k)] + c_moments[I3D(i-1,j+1,k)] - 
                 c_moments[I3D(i,j-1,k)] - c_moments[I3D(i-1,j-1,k)]) / (4.0 * dy);
            
            double dx_ujp = 1.0 / sigma_jp * 
                (c_moments[I3D(i+1,j+1,k)] + c_moments[I3D(i+1,j,k)] - 
                 c_moments[I3D(i-1,j+1,k)] - c_moments[I3D(i-1,j,k)]) / (4.0 * dx);
            double dx_ujn = 1.0 / sigma_jn * 
                (c_moments[I3D(i+1,j,k)] + c_moments[I3D(i+1,j-1,k)] - 
                 c_moments[I3D(i-1,j,k)] - c_moments[I3D(i-1,j-1,k)]) / (4.0 * dx);
            double dy_ujp = 1.0 / sigma_jp * 
                (c_moments[I3D(i,j+1,k)] - c_moments[I3D(i,j,k)]) / dy;
            double dy_ujn = 1.0 / sigma_jn * 
                (c_moments[I3D(i,j,k)] - c_moments[I3D(i,j-1,k)]) / dy;
            
            for(int k2 = 0; k2 < c_numMoments; k2++)
            {
                flux[I3D(i,j,k2)] += (dx_uip - dx_uin) / dx * 
                    (c_flux2OperatorXiXi[k2*c_numMoments + k] - 
                     c_fluxOperatorXiXi[k2*c_numMoments + k]);
                flux[I3D(i,j,k2)] += (dy_uip - dy_uin) / dx * 
                    (c_flux2OperatorXiEta[k2*c_numMoments + k] - 
                     c_fluxOperatorXiEta[k2*c_numMoments + k]);
                flux[I3D(i,j,k2)] += (dx_ujp - dx_ujn) / dy * 
                    (c_flux2OperatorEtaXi[k2*c_numMoments + k] - 
                     c_fluxOperatorXiEta[k2*c_numMoments + k]);
                flux[I3D(i,j,k2)] += (dy_ujp - dy_ujn) / dy * 
                    (c_flux2OperatorEtaEta[k2*c_numMoments + k] - 
                     c_fluxOperatorEtaEta[k2*c_numMoments + k]);
            }
        }
    }

    #ifdef USE_PAPI
    profile_finish_update("DnSolver::solveFlux");
    #endif
}


/*
    Updates one time step.
*/
void DnSolver::update(double dt, double dx, double dy)
{
    #ifdef USE_PAPI
    profile_start_update("DnSolver::update");
    #endif

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
                if(c_initCond == 1)
                {
                    double x_i = c_initX + (i-c_gX[1]) * dx;
                    double y_j = c_initY + (j-c_gY[1]) * dy;
                    if(fabs(x_i) < 0.5 && fabs(y_j) < 0.5)
                        c_moments[I3D(i,j,k)] = c_moments[I3D(i,j,k)] + dt * 2.0 * sqrt(M_PI);
                }
            }
        }
    }


    // Do first Euler step.
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
                if(c_initCond == 1)
                {
                    double x_i = c_initX + (i-c_gX[1]) * dx;
                    double y_j = c_initY + (j-c_gY[1]) * dy;
                    if(fabs(x_i) < 0.5 && fabs(y_j) < 0.5)
                        c_moments[I3D(i,j,k)] = c_moments[I3D(i,j,k)] + dt * 2.0 * sqrt(M_PI);
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

    #ifdef USE_PAPI
    profile_finish_update("DnSolver::update");
    #endif
}

