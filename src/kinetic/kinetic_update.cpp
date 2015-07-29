/*
 File:   update.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "kinetic_solver.h"
#include "../utils.h"
#include <strings.h>
#include <math.h>
#include <stdlib.h>

#ifdef USE_PAPI
#include "../profiling.h"
#endif

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
void KineticSolver::solveFlux(double *kinetic, double *flux, double dx, double dy)
{
    #ifdef USE_PAPI
    profile_start_update("KineticSolver::solveFlux");
    #endif
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

        
        // Add up each quadrature point's flux.
        for(int q = 0; q < c_numQuadPoints; q++)
        {
            double northFlux;
            double southFlux;
            double eastFlux;
            double westFlux;
            
            if(c_eta[q] > 0)
            {
                double s1 = kinetic[I3D(i,j-2,q)];
                double s2 = kinetic[I3D(i,j-1,q)];
                double s3 = kinetic[I3D(i,j,q)];
                double s4 = kinetic[I3D(i,j+1,q)];
                northFlux = s3 + 0.5 * slopefit(s2, s3, s4, 2.0);
                southFlux = s2 + 0.5 * slopefit(s1, s2, s3, 2.0);
            }
            else
            {
                double s1 = kinetic[I3D(i,j-1,q)];
                double s2 = kinetic[I3D(i,j,q)];
                double s3 = kinetic[I3D(i,j+1,q)];
                double s4 = kinetic[I3D(i,j+2,q)];
                northFlux = s3 - 0.5 * slopefit(s2, s3, s4, 2.0);
                southFlux = s2 - 0.5 * slopefit(s1, s2, s3, 2.0);
            }
            if(c_xi[q] > 0)
            {
                double s1 = kinetic[I3D(i-2,j,q)];
                double s2 = kinetic[I3D(i-1,j,q)];
                double s3 = kinetic[I3D(i,j,q)];
                double s4 = kinetic[I3D(i+1,j,q)];
                eastFlux = s3 + 0.5 * slopefit(s2, s3, s4, 2.0);
                westFlux = s2 + 0.5 * slopefit(s1, s2, s3, 2.0);
            }
            else
            {
                double s1 = kinetic[I3D(i-1,j,q)];
                double s2 = kinetic[I3D(i,j,q)];
                double s3 = kinetic[I3D(i+1,j,q)];
                double s4 = kinetic[I3D(i+2,j,q)];
                eastFlux = s3 - 0.5 * slopefit(s2, s3, s4, 2.0);
                westFlux = s2 - 0.5 * slopefit(s1, s2, s3, 2.0);
            }
            
            flux[I3D(i,j,q)] = (eastFlux - westFlux) * c_xi[q] / dx + 
                               (northFlux - southFlux) * c_eta[q] / dy;
        }
    }
    #ifdef USE_PAPI
    profile_finish_update("KineticSolver::solveFlux");
    #endif
}


/*
    Updates one time step.
*/
void KineticSolver::update(double dt, double dx, double dy)
{
    static bool firstTime = true;
    static double *kineticOld;
    if(firstTime)
    {
        firstTime = false;
        kineticOld = (double*)malloc((c_gX[3]-c_gX[0]+1) * (c_gY[3]-c_gY[0]+1) * 
            c_numQuadPoints * sizeof(double));
    }
    
    #ifdef USE_PAPI
    profile_start_update("KineticSolver::update");
    #endif
    
    // Copy initial grid for use later.
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                kineticOld[I3D(i,j,q)] = c_kinetic[I3D(i,j,q)];
            }
        }
    }
    

    // Do first Euler step.
    communicateBoundaries();
    solveFlux(c_kinetic, c_flux, dx, dy);

    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            // Calculate integral of F.
            double integral = 0;
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                double f = c_kinetic[I3D(i,j,q)];
                double w =  c_quadWeights[q];
                integral += w * f;
            }
            integral = integral * c_sigmaS[I2D(i,j)] / (4 * M_PI);

            // Do Euler Step.
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                c_kinetic[I3D(i,j,q)] = c_kinetic[I3D(i,j,q)] - 
                    dt * c_sigmaT[I2D(i,j)] * c_kinetic[I3D(i,j,q)] + dt * integral;
                if(c_initCond == 1)
                {
                    double x_i = c_initX + i * dx;
                    double y_j = c_initY + j * dy;
                    if(fabs(x_i) < 0.5 && fabs(y_j) < 0.5)
                        c_kinetic[I3D(i,j,q)] = c_kinetic[I3D(i,j,q)] + dt;
                }
                c_kinetic[I3D(i,j,q)] = c_kinetic[I3D(i,j,q)] - dt * c_flux[I3D(i,j,q)];
            }
        }
    }
    
    
    // Do second Euler step.
    communicateBoundaries();
    solveFlux(c_kinetic, c_flux, dx, dy);

    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            // Calculate integral of F.
            double integral = 0;
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                double f = c_kinetic[I3D(i,j,q)];
                double w =  c_quadWeights[q];
                integral += w * f;
            }
            integral = integral * c_sigmaS[I2D(i,j)] / (4 * M_PI);

            // Do Euler Step.
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                c_kinetic[I3D(i,j,q)] = c_kinetic[I3D(i,j,q)] - 
                    dt * c_sigmaT[I2D(i,j)] * c_kinetic[I3D(i,j,q)] + dt * integral;
                if(c_initCond == 1)
                {
                    double x_i = c_initX + i * dx;
                    double y_j = c_initY + j * dy;
                    if(fabs(x_i) < 0.5 && fabs(y_j) < 0.5)
                        c_kinetic[I3D(i,j,q)] = c_kinetic[I3D(i,j,q)] + dt;
                }
                c_kinetic[I3D(i,j,q)] = c_kinetic[I3D(i,j,q)] - dt * c_flux[I3D(i,j,q)];
            }
        }
    }
    
    
    // Average initial with second Euler step.
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                c_kinetic[I3D(i,j,q)] = 0.5 * (c_kinetic[I3D(i,j,q)] + kineticOld[I3D(i,j,q)]);
            }
        }
    }

    #ifdef USE_PAPI
    profile_finish_update("KineticSolver::update");
    #endif
}

