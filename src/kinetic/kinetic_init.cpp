/*
 File:   init.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "kinetic_solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"


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
double KineticSolver::init(double dx, double dy)
{
    // Read moment.deck
    FILE *file = fopen("kinetic.deck", "r");
    if(file == NULL)
    {
        printf("Could not open input file: %s\n", "moment.deck");
        utils_abort();
    }

    int quadOrder;
    double cflFactor;
    utils_readLine(file, &quadOrder);
    utils_readLine(file, &cflFactor);

    fclose(file);


    // Print out moment.deck
    if(c_node == 0)
    {
        printf("quadOrder:   %d\n", quadOrder);
        printf("cflFactor:   %f\n", cflFactor);
    }
    
    
    // Maximum value for delta t.
    double maxDt = cflFactor * (2.0 / (2.0 + 2.0)) * (dx * dy) / (dx + dy);


    // Number of quadrature points.
    c_numQuadPoints = quadOrder * quadOrder;
    c_vectorSize    = c_numQuadPoints;
    
    
    // The weights and nodes for integration and the angles.
    double *w   = (double*)malloc(quadOrder * sizeof(double));
    double *mu  = (double*)malloc(quadOrder * sizeof(double));
    double *phi = (double*)malloc(2 * quadOrder * sizeof(double));
    
    
    // Allocate memory for quadrature.
    c_quadWeights = (double*)malloc(c_numQuadPoints * sizeof(double));
    c_xi = (double*)malloc(c_numQuadPoints * sizeof(double));
    c_eta = (double*)malloc(c_numQuadPoints * sizeof(double));
    
    
    // Allocate memory for grid cells.
    int numXCells = 2 * NUM_GHOST_CELLS + c_sizeX;
    int numYCells = 2 * NUM_GHOST_CELLS + c_sizeY;
    c_kinetic = (double*)malloc(numXCells * numYCells * c_numQuadPoints * sizeof(double));
    c_flux    = (double*)malloc(numXCells * numYCells * c_numQuadPoints * sizeof(double));
    
    
    // Initial conditions
    for(int i = c_gX[0]; i <= c_gX[3]; i++)
    {
        for(int j = c_gY[0]; j <= c_gY[3]; j++)
        {
            for(int q = 0; q < c_numQuadPoints; q++)
            {
                c_kinetic[I3D(i,j,q)] = c_initialGrid[I2D(i,j)];
            }
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
            double muTemp = mu[q1];
            double phiTemp = phi[q2];
            double xiTemp = sqrt(1.0 - muTemp * muTemp) * cos(phiTemp);
            double etaTemp = sqrt(1.0 - muTemp * muTemp) * sin(phiTemp);
            int q = 2 * q1 * quadOrder + q2;     // quadrature counter
            c_quadWeights[q] = wTemp;
            c_xi[q] = xiTemp;
            c_eta[q] = etaTemp;
        }
    }


    // Free memory.
    free(w);
    free(mu);
    free(phi);

    return maxDt;
}

