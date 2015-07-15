/*
 File:   utils.cpp
 Author: Kris Garrett
 Date:   July 26, 2012
 
 Several functions are defined here as helper utilities to the rest of the 
 program.

*/


#include "utils.h"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "lib/sphere_lebedev_rule.h"

#ifdef USE_MPI
#include <mpi.h>
#endif


/*
 Implements abort for both serial and MPI code.
*/
extern "C"
void utils_abort()
{
    #ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    #else
    abort();
    #endif
}


double utils_norm1(int n, double *A)
{
    double oneNorm = 0.0;
    
    for(int i = 0; i < n; i++)
    {
        double oneNormTemp = 0;
        for(int j = 0; j < n; j++)
        {
            oneNormTemp += fabs(A[i*n+j]);
        }
        if(oneNormTemp > oneNorm)
            oneNorm = oneNormTemp;
    }
    
    return oneNorm;
}


void utils_getGaussianWeightsAndNodes(int n, double *w, double *mu)
{
    gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(n);
    for(int i = 0; i < n; i++)
        gsl_integration_glfixed_point(-1, 1, i, &mu[i], &w[i], table);
}

void utils_getLebedevWeightsAndNodes(int numPoints, double *w, double *x, double *y)
{
    double *z = (double *)malloc(numPoints * sizeof(double));
    ld_by_order(numPoints, x, y, z, w);
    free(z);

    for(int i = 0; i < numPoints; i++) {
        w[i] *= 4 * M_PI;
    }
}

int utils_numLebedevQuadPoints(int rule_number) {
    if(rule_number < 1 || rule_number > 65) {
        printf("Choose 0 < QUAD_ORDER < 66\n");
        utils_abort();
    }
    if(available_table(rule_number) == 1) {
        return order_table(rule_number);
    } else {
        printf("Invalid QUAD_ORDER\n");
        utils_abort();
        return -1;
    }
}

void test_quadature(int numPoints, double *w, double *x, double *y) {
    printf("Testing quadrature\n");

    double integral = 0.0;
    for(int i = 0; i < numPoints; i++) {
        integral += w[i];
    }
    printf("<1> = %f\t\terror %f\n", integral, fabs(integral - 4 * M_PI));

    integral = 0.0;
    for(int i = 0; i < numPoints; i++) {
        integral += w[i] * x[i];
    }
    printf("<x> = %f\t\terror %f\n", integral, fabs(integral));

    integral = 0.0;
    for(int i = 0; i < numPoints; i++) {
        integral += w[i] * x[i] * x[i];
    }
    printf("<x^2> = %f\terror %f\n", integral, fabs(integral - 4 * M_PI / 3));
}

