/*
 File:   init.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/


#include "momopt_solver.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../input_deck_reader.h"
#include "opt/opt.h"
#include "opt/fobj.h"
#include <gsl/gsl_linalg.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

/*
    Checks the status of input parameters.
*/
static void checkInput(bool isOk, int lineNumber)
{
    if(!isOk) {
        printf("momopt_init.cpp:input error at line %d\n", lineNumber);
        utils_abort();
    }
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
double MomOptSolver::init(double dx, double dy)
{
    int momentOrder, quadOrder;
    double cflFactor;

    checkInput(c_inputDeckReader.getValue("MOMOPT_MOMENT_ORDER", &momentOrder), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_QUAD_ORDER", &quadOrder), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_CFL_FACTOR", &cflFactor), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_TOL", &c_tolerance), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_COND_H_MAX", &c_condHMax), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_COND_H_MAX_BFGS", &c_condHMaxBfgs), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_MAX_ITER", &c_maxIter), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_MAX_BFGS_ITER", &c_maxIterBfgs), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_USE_GAUNT", &c_useGaunt), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_THETA", &c_theta), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_DELTA_PPN", &c_deltaPPn), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_MOMENT_TYPE", &c_momentType), __LINE__);
    checkInput(c_inputDeckReader.getValue("MOMOPT_OPTIMIZATION_TYPE", &c_optType), __LINE__);
    
    
    // Maximum value for delta t.
    double maxDt = cflFactor * (2.0 / (2.0 + c_theta)) * (dx * dy) / (dx + dy);


    // Optimization Statistics
    c_optStats.maxIter = 0;
    //c_optStats.maxGammaIter = 0;
    c_optStats.iterMean = 0.0;
    //c_optStats.iterGammaMean = 0.0;
    c_optStats.numDualSolves = 0;
    for(int i = 0; i < NUM_REGULARIZATIONS; i++)
    {
        c_optStats.histReg[i] = 0;
    }
    
    
    // Number of quadrature points and number of moments.
    c_numQuadPoints = quadOrder * quadOrder;
    c_numMoments    = (momentOrder + 1) * (momentOrder + 2) / 2;
    c_numManyMoments = (momentOrder + 1) * (momentOrder + 2) / 2;
    c_vectorSize    = c_numMoments;
    if(c_useGaunt == 1)
        c_numManyMoments = (2 * momentOrder + 1) * (2 * momentOrder + 2) / 2;
    
    
    // The weights and nodes for integration and the angles.
    gsl_vector *w   = gsl_vector_calloc(quadOrder);
    gsl_vector *mu  = gsl_vector_calloc(quadOrder);
    gsl_vector *phi = gsl_vector_calloc(2 * quadOrder);
    
    
    // Allocate memory for quadrature.
    c_spHarm = (double*)malloc(c_numQuadPoints * c_numManyMoments * sizeof(double));
    c_quadWeights = (double*)malloc(c_numQuadPoints * sizeof(double));
    c_xi  = (double*)malloc(c_numQuadPoints * sizeof(double));
    c_eta = (double*)malloc(c_numQuadPoints * sizeof(double));
    
    
    // Allocate memory for grid cells.
    int numXCells = 2 * NUM_GHOST_CELLS + c_sizeX;
    int numYCells = 2 * NUM_GHOST_CELLS + c_sizeY;
    c_moments = (double*)malloc(numXCells * numYCells * c_numMoments * sizeof(double));
    c_flux    = (double*)malloc(numXCells * numYCells * c_numMoments * sizeof(double));
    c_alpha   = (double*)malloc(numXCells * numYCells * c_numMoments * sizeof(double));
    c_alphaP  = (double*)malloc(numXCells * numYCells * c_numMoments * sizeof(double));
    c_P2M     = (double*)malloc(numXCells * numYCells * c_numMoments * c_numMoments * sizeof(double));
    
    
    // Initial conditions
    for(int i = c_gX[0]; i <= c_gX[3]; i++)
    {
        for(int j = c_gY[0]; j <= c_gY[3]; j++)
        {
            // Set moments
            c_moments[I3D(i,j,0)] = c_initialGrid[I2D(i,j)] * 2.0 * sqrt(M_PI);
            for(int k = 1; k < c_numMoments; k++)
                c_moments[I3D(i,j,k)] = 0.0;
            
            
            // Set alpha
            for(int k = 1; k < c_numMoments; k++)
                c_alpha[I3D(i,j,k)] = 0.0;
            if(c_momentType == MOMENT_TYPE_MN)
            {
                c_alpha[I3D(i,j,0)] = 2.0 * sqrt(M_PI) * 
                    log(c_moments[I3D(i,j,0)] / (2.0 * sqrt(M_PI)));
            }
            else if(c_momentType == MOMENT_TYPE_PPN)
            {
                c_alpha[I3D(i,j,0)] = c_moments[I3D(i,j,0)] - 
                    4.0 * c_deltaPPn * M_PI / c_moments[I3D(i,j,0)];
            }
            
            
            // Set Change of Basis stuff.
            if(c_optType == OPT_TYPE_CHANGE_BASIS || c_optType == OPT_TYPE_BFGS)
            {
                for(int k = 1; k < c_numMoments; k++)
                    c_alphaP[I3D(i,j,k)] = 0.0;
                if(c_momentType == MOMENT_TYPE_MN)
                    c_alphaP[I3D(i,j,0)] = 2.0 * sqrt(M_PI) * log(0.5 / sqrt(M_PI));
                else if(c_momentType == MOMENT_TYPE_PPN)
                    c_alphaP[I3D(i,j,0)] = 1.0 - 4.0 * c_deltaPPn * M_PI;
                
                memset(&c_P2M[I3D2(i,j,0)], 0, c_numMoments * c_numMoments * sizeof(double));
                for(int k = 0; k < c_numMoments; k++)
                    c_P2M[I3D2(i,j,k)] = 1.0;
            }
        }
    }
    
    
    // Get quadrature.
    utils_getGaussianWeightsAndNodes(quadOrder, w->data, mu->data);
    
    
    // Azimuthal angles of quadrature.
    for(int k = 0; k < 2 * quadOrder; k++)
    {
        gsl_vector_set(phi, k, (k + 0.5) * M_PI / quadOrder);
    }
    
    
    // Setup quadrature.
    // Gaussian quadrature on z-axis.
    for(int q1 = 0; q1 < quadOrder / 2; q1++)
    {
        double wTemp = 2.0 * M_PI / quadOrder * gsl_vector_get(w, q1);
        for(int q2 = 0; q2 < 2 * quadOrder; q2++)
        {
            double shTemp;
            double muTemp = gsl_vector_get(mu, q1);
            double phiTemp = gsl_vector_get(phi, q2);
            double xiTemp = sqrt(1.0 - muTemp * muTemp) * cos(phiTemp);
            double etaTemp = sqrt(1.0 - muTemp * muTemp) * sin(phiTemp);
            int k = 0;                          // moment counter
            int q = 2 * q1 * quadOrder + q2;     // quadrature counter
            c_quadWeights[q] = wTemp;
            c_xi[q] = xiTemp;
            c_eta[q] = etaTemp;

            // distribute measure at each level equally among 2N azimuthal points
            int N = momentOrder + 1;
            if(c_useGaunt == 1)
                N = 2 * momentOrder + 1;
            for(int n = 0; n < N; n++)
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

                        c_spHarm[q*c_numManyMoments+k] = shTemp;
                        k++;
                    }
                }
            }
        }
    }
    
    
    // Solve for the Gaunt coefficients.
    if(c_useGaunt == 1)
    {
        g_hStruct = new HStruct[c_numMoments * c_numMoments];
        for(int i = 0; i < c_numMoments * c_numMoments; i++)
            g_hStruct[i].s = new HStructData[5 * c_numMoments];

        gsl_vector *b = gsl_vector_alloc(c_numQuadPoints);
        gsl_vector *beta = gsl_vector_alloc(c_numManyMoments);
        gsl_vector *tau = gsl_vector_alloc(c_numManyMoments);
        gsl_vector *residual = gsl_vector_alloc(c_numQuadPoints);
        gsl_matrix *A = gsl_matrix_alloc(c_numQuadPoints, c_numManyMoments);
        memcpy(A->data, c_spHarm, c_numQuadPoints * c_numManyMoments * sizeof(double));
        gsl_linalg_QR_decomp(A, tau);
        
        for(int i = 0; i < c_numMoments; i++)
        {
            for(int j = 0; j < c_numMoments; j++)
            {
                for(int q = 0; q < c_numQuadPoints; q++)
                {
                    gsl_vector_set(b, q, c_spHarm[q*c_numManyMoments+i] * c_spHarm[q*c_numManyMoments+j]);
                }
                gsl_linalg_QR_lssolve(A, tau, b, beta, residual);
                
                int index = i * c_numMoments + j;
                g_hStruct[index].numMoments = 0;
                for(int k = 0; k < c_numManyMoments; k++)
                {
                    if(fabs(gsl_vector_get(beta, k)) > 10e-10)
                    {
                        int momentIndex = g_hStruct[index].numMoments;
                        g_hStruct[index].numMoments++;
                        g_hStruct[index].s[momentIndex].moment = k;
                        g_hStruct[index].s[momentIndex].value = gsl_vector_get(beta, k);
                    }
                }
            }
        }
    }


    // Free memory.
    gsl_vector_free(w);
    gsl_vector_free(mu);
    gsl_vector_free(phi);

    // Let functions initialize.
    initUpdate(c_numMoments);
    initFobj(c_numOmpThreads, c_numManyMoments);
    
    return maxDt;
}

