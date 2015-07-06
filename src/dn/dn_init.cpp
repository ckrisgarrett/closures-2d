/*
 File:   init.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "dn_solver.h"
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../input_deck_reader.h"

/*
    Checks the status of input parameters.
*/
static void checkInput(bool isOk, int lineNumber)
{
    if(!isOk) {
        printf("dn_init.cpp:input error at line %d\n", lineNumber);
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
double DnSolver::init(double dx, double dy)
{
    int momentOrder, quadOrder;
    double cflFactor;

    checkInput(c_inputDeckReader.getValue("MOMENT_ORDER", &momentOrder), __LINE__);
    checkInput(c_inputDeckReader.getValue("QUAD_ORDER", &quadOrder), __LINE__);
    checkInput(c_inputDeckReader.getValue("CFL_FACTOR", &cflFactor), __LINE__);
    
    
    // Maximum value for delta t.
    double maxDt = cflFactor * 3.0 / 4.0 * MIN(dx, dy) * MIN(dx, dy);


    // Number of quadrature points and number of moments.
    c_numQuadPoints = quadOrder * quadOrder;
    c_numMoments    = (momentOrder + 1) * (momentOrder + 2) / 2;
    c_vectorSize    = c_numMoments;
    
    
    // The weights and nodes for integration and the angles.
    double *w   = (double*)malloc(quadOrder * sizeof(double));
    double *mu  = (double*)malloc(quadOrder * sizeof(double));
    double *phi = (double*)malloc(2 * quadOrder * sizeof(double));
    
    
    // Allocate memory for quadrature.
    c_p = (double*)malloc(c_numQuadPoints * c_numMoments * sizeof(double));
    c_w = (double*)malloc(c_numQuadPoints * sizeof(double));
    c_xi  = (double*)malloc(c_numQuadPoints * sizeof(double));
    c_eta = (double*)malloc(c_numQuadPoints * sizeof(double));
    
    
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
            int q = 2 * q1 * quadOrder + q2;     // quadrature counter
            c_w[q] = wTemp;
            c_xi[q] = xiTemp;
            c_eta[q] = etaTemp;

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

                        c_p[q*c_numMoments+k] = shTemp;
                        k++;
                    }
                }
            }
        }
    }
    
    
    // Set up matrices for flux operators.
    c_fluxOperatorXi       = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_fluxOperatorEta      = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_fluxOperatorXiPlus   = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_fluxOperatorXiMinus  = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_fluxOperatorEtaPlus  = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_fluxOperatorEtaMinus = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_fluxOperatorXiXi     = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_fluxOperatorXiEta    = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_fluxOperatorEtaEta   = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    for(int i = 0; i < c_numMoments; i++)
    {
        for(int j = 0; j < c_numMoments; j++)
        {
            c_fluxOperatorXi[i*c_numMoments+j]       = 0.0;
            c_fluxOperatorEta[i*c_numMoments+j]      = 0.0;
            c_fluxOperatorXiPlus[i*c_numMoments+j]   = 0.0;
            c_fluxOperatorXiMinus[i*c_numMoments+j]  = 0.0;
            c_fluxOperatorEtaPlus[i*c_numMoments+j]  = 0.0;
            c_fluxOperatorEtaMinus[i*c_numMoments+j] = 0.0;
            c_fluxOperatorXiXi[i*c_numMoments+j]     = 0.0;
            c_fluxOperatorXiEta[i*c_numMoments+j]    = 0.0;
            c_fluxOperatorEtaEta[i*c_numMoments+j]   = 0.0;
            
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                double m_i   = c_p[q*c_numMoments+i];
                double m_j   = c_p[q*c_numMoments+j];
                double w_q   = c_w[q];
                
                c_fluxOperatorXi[i*c_numMoments+j]     += c_xi[q] * m_i * m_j * w_q;
                c_fluxOperatorEta[i*c_numMoments+j]    += c_eta[q] * m_i * m_j * w_q;
                c_fluxOperatorXiXi[i*c_numMoments+j]   += c_xi[q] * c_xi[q]  * m_i * m_j * w_q;
                c_fluxOperatorXiEta[i*c_numMoments+j]  += c_xi[q] * c_eta[q] * m_i * m_j * w_q;
                c_fluxOperatorEtaEta[i*c_numMoments+j] += c_eta[q] * c_eta[q] * m_i * m_j * w_q;
                
                if(c_xi[q] > 0)
                    c_fluxOperatorXiPlus[i*c_numMoments+j] += c_xi[q] * m_i * m_j * w_q;
                else
                    c_fluxOperatorXiMinus[i*c_numMoments+j] += c_xi[q] * m_i * m_j * w_q;
                if(c_eta[q] > 0)
                    c_fluxOperatorEtaPlus[i*c_numMoments+j] += c_eta[q] * m_i * m_j * w_q;
                else
                    c_fluxOperatorEtaMinus[i*c_numMoments+j] += c_eta[q] * m_i * m_j * w_q;
            }
        }
    }
    
    
    char nChar = 'N';
    double one = 1.0;
    double zero = 0.0;
    c_flux2OperatorXiXi   = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_flux2OperatorXiEta  = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_flux2OperatorEtaXi  = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    c_flux2OperatorEtaEta = (double*)malloc(c_numMoments * c_numMoments * sizeof(double));
    dgemm_(&nChar, &nChar, &c_numMoments, &c_numMoments, &c_numMoments, &one, 
        c_fluxOperatorXi, &c_numMoments, c_fluxOperatorXi, &c_numMoments, 
        &zero, c_flux2OperatorXiXi, &c_numMoments);
    dgemm_(&nChar, &nChar, &c_numMoments, &c_numMoments, &c_numMoments, &one, 
        c_fluxOperatorXi, &c_numMoments, c_fluxOperatorEta, &c_numMoments, 
        &zero, c_flux2OperatorXiEta, &c_numMoments);
    dgemm_(&nChar, &nChar, &c_numMoments, &c_numMoments, &c_numMoments, &one, 
        c_fluxOperatorEta, &c_numMoments, c_fluxOperatorXi, &c_numMoments, 
        &zero, c_flux2OperatorEtaXi, &c_numMoments);
    dgemm_(&nChar, &nChar, &c_numMoments, &c_numMoments, &c_numMoments, &one, 
        c_fluxOperatorEta, &c_numMoments, c_fluxOperatorEta, &c_numMoments, 
        &zero, c_flux2OperatorEtaEta, &c_numMoments);
    
    
    // Free memory.
    free(w);
    free(mu);
    free(phi);


    return maxDt;
}

