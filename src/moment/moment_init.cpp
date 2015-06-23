/*
 File:   init.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "moment_solver.h"
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"


/*
    The filter used for the FPn2 method as described by Radice et. al.
*/
static double filterSSpline(double x)
{
    return 1.0 / (1.0 + x * x * x * x);
}

static double filterLanczos(double x)
{
    if(x < .0000001)
        return 1.0;
    return sin(x) / x;
}


/*
    This function initializes the data.
    
    node:       MPI node index (0 if not using MPI)
    numNodes:   number of MPI nodes (1 if not using MPI)
    dx:         delta x
    dy:         delta y
    
    Returns maximum for dt.
    
    Notes:
    Must be run after setting up MPI and OpenMP if using those options.
*/
double MomentSolver::init(double dx, double dy)
{
    // Read moment.deck
    FILE *file = fopen("moment.deck", "r");
    if(file == NULL)
    {
        printf("Could not open input file: %s\n", "moment.deck");
        utils_abort();
    }

    int momentOrder, quadOrder;
    double cflFactor;
    utils_readLine(file, &momentOrder);
    utils_readLine(file, &quadOrder);
    utils_readLine(file, &c_filterType);
    utils_readLine(file, &c_filterTune);
    utils_readLine(file, &cflFactor);

    fclose(file);


    // Print out moment.deck
    if(c_node == 0)
    {
        printf("momentOrder: %d\n", momentOrder);
        printf("quadOrder:   %d\n", quadOrder);
        printf("filterType:  %d\n", c_filterType);
        printf("filterTune:  %f\n", c_filterTune);
        printf("cflFactor:   %f\n", cflFactor);
    }
    
    
    // Maximum value for delta t.
    double maxDt = cflFactor * 0.5 * MIN(dx, dy);


    // Number of quadrature points and number of moments.
    c_numQuadPoints = quadOrder * quadOrder;
    c_numMoments    = (momentOrder + 1) * (momentOrder + 2) / 2;
    c_vectorSize    = c_numMoments;
    
    
    // The weights and nodes for integration and the angles.
    double *w   = (double*)malloc(quadOrder * sizeof(double));
    double *mu  = (double*)malloc(quadOrder * sizeof(double));
    double *phi = (double*)malloc(2 * quadOrder * sizeof(double));
    
    
    // Allocate memory for quadrature.
    double *spHarm = (double*)malloc(c_numQuadPoints * c_numMoments * sizeof(double));
    double *quadWeights = (double*)malloc(c_numQuadPoints * sizeof(double));
    double *xi  = (double*)malloc(c_numQuadPoints * sizeof(double));
    double *eta = (double*)malloc(c_numQuadPoints * sizeof(double));
    
    
    // Allocate memory for flux operators
    c_pnFluxOperatorXiPlus   = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_pnFluxOperatorXiMinus  = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_pnFluxOperatorEtaPlus  = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_pnFluxOperatorEtaMinus = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    
    
    // Allocate memory for grid cells.
    int numXCells = 2 * NUM_GHOST_CELLS + c_sizeX;
    int numYCells = 2 * NUM_GHOST_CELLS + c_sizeY;
    c_moments = (double*)malloc(numXCells * numYCells * c_numMoments * sizeof(double));
    c_flux    = (double*)malloc(numXCells * numYCells * c_numMoments * sizeof(double));
    
    
    // Initial conditions
    for(int i = c_gX[0]; i <= c_gX[3]; i++)
    {
        for(int j = c_gY[0]; j <= c_gY[3]; j++)
        {
            c_moments[I3D(i,j,0)] = c_initialGrid[I2D(i,j)] * 2.0 * sqrt(M_PI);
            for(int k = 1; k < c_numMoments; k++)
                c_moments[I3D(i,j,k)] = 0.0;
        }
    }
    
    
    // Get quadrature.
    utils_getGaussianWeightsAndNodes(quadOrder, w, mu);
    
    
    // Azimuthal angles of quadrature.
    for(int k = 0; k < 2 * quadOrder; k++)
    {
        phi[k] = (k + 0.5) * M_PI / quadOrder;
    }
    
    
    // Setup quadrature.
    // Gaussian quadrature on z-axis.
    for(int q1 = 0; q1 < quadOrder / 2; q1++)
    {
        double wTemp = 2.0 * M_PI / quadOrder * w[q1];
        for(int q2 = 0; q2 < 2 * quadOrder; q2++)
        {
            double shTemp;
            double muTemp = mu[q1];
            double phiTemp = phi[q2];
            double xiTemp = sqrt(1.0 - muTemp * muTemp) * cos(phiTemp);
            double etaTemp = sqrt(1.0 - muTemp * muTemp) * sin(phiTemp);
            int k = 0;                          // moment counter
            int q = 2 * q1 * quadOrder + q2;    // quadrature counter
            quadWeights[q] = wTemp;
            xi[q] = xiTemp;
            eta[q] = etaTemp;
            
            // distribute measure at each level equally among 2N azimuthal points
            for(int n = 0; n < momentOrder + 1; n++)
            {
                for(int m = -n; m <= n; m++)
                {
                    if((m + n) % 2 == 0)
                    {
                        if(m < 0)
                            shTemp = M_SQRT2 * sin(m * phiTemp) * 
                                     gsl_sf_legendre_sphPlm(n, -m, muTemp);
                        else if(m > 0)
                            shTemp = M_SQRT2 * cos(m * phiTemp) * 
                                     gsl_sf_legendre_sphPlm(n, m, muTemp);
                        else
                            shTemp = gsl_sf_legendre_sphPlm(n, 0, muTemp);

                        spHarm[q*c_numMoments+k] = shTemp;
                        k++;
                    }
                }
            }
        }
    }
    
    
    // Set up matrices for Pn flux operators.
    for(int i = 0; i < c_numMoments; i++)
    {
        for(int j = 0; j < c_numMoments; j++)
        {
            c_pnFluxOperatorXiPlus[i*c_numMoments+j] = 0.0;
            c_pnFluxOperatorXiMinus[i*c_numMoments+j] = 0.0;
            c_pnFluxOperatorEtaPlus[i*c_numMoments+j] = 0.0;
            c_pnFluxOperatorEtaMinus[i*c_numMoments+j] = 0.0;
            
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                double m_i   = spHarm[q*c_numMoments+i];
                double m_j   = spHarm[q*c_numMoments+j];
                double w_q   = quadWeights[q];
                
                double valueXi  = xi[q] * m_i * m_j * w_q;
                double valueEta = eta[q] * m_i * m_j * w_q;
                
                if(xi[q] > 0)
                {
                    c_pnFluxOperatorXiPlus[i*c_numMoments+j] += valueXi;
                }
                else
                {
                    c_pnFluxOperatorXiMinus[i*c_numMoments+j] += valueXi;
                }
                if(eta[q] > 0)
                {
                    c_pnFluxOperatorEtaPlus[i*c_numMoments+j] += valueEta;
                }
                else
                {
                    c_pnFluxOperatorEtaMinus[i*c_numMoments+j] += valueEta;
                }
            }
        }
    }

    
    // Create g_momentFilter.
    c_momentFilter = (double*)malloc(c_numMoments * sizeof(double));
    switch(c_filterType)
    {
        case FILTER_TYPE_NONE:
        {
            for(int k = 0; k < c_numMoments; k++)
            {
                c_momentFilter[k] = 1.0;
            }
            break;
        }
        case FILTER_TYPE_HAUCK:
        {
            int nm = 0;
            for(int n = 0; n <= momentOrder; n++)
            {
                for(int m = -n; m <= n; m++)
                {
                    if((m + n) % 2 == 0)
                    {
                        double omega = c_filterTune;
                        double N = momentOrder;
                        double L = 1.0;
                        double sigma_t = 1.0;
                        double alpha = omega / (N * N) / (sigma_t * L + N) / (sigma_t * L + N);
                        c_momentFilter[nm] = 1.0 / (1.0 + alpha * n * n * (n+1) * (n+1));
                        nm++;
                    }
                }
            }
            break;
        }
        case FILTER_TYPE_SSPLINE:
        {
            int nm = 0;
            for(int n = 0; n <= momentOrder; n++)
            {
                for(int m = -n; m <= n; m++)
                {
                    if((m + n) % 2 == 0)
                    {
                        double N = momentOrder;
                        double sigma_e = c_filterTune;
                        double s = -sigma_e * maxDt / 
                            log(filterSSpline(N / (N+1)));
                        c_momentFilter[nm] = pow(filterSSpline(n / (N + 1)), s);
                        nm++;
                    }
                }
            }
            break;
        }
        case FILTER_TYPE_LANCZOS:
        {
            int nm = 0;
            for(int n = 0; n <= momentOrder; n++)
            {
                for(int m = -n; m <= n; m++)
                {
                    if((m + n) % 2 == 0)
                    {
                        double N = momentOrder;
                        double sigma_e = c_filterTune;
                        double s = -sigma_e * maxDt / 
                            log(filterLanczos(N / (N+1)));
                        c_momentFilter[nm] = pow(filterLanczos(n / (N + 1)), s);
                        nm++;
                    }
                }
            }
            break;
        }
        default:
            break;
    }


    // Free memory.
    free(w);
    free(mu);
    free(phi);
    free(spHarm);
    free(quadWeights);
    free(xi);
    free(eta);


    return maxDt;
}

