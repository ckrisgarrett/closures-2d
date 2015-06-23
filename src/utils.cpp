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


/*
 Reads a line of a file of the form <Name> <Value>
*/
void utils_readLine(FILE *file, double *x)
{
    int readReturn;
    char uselessString[1024];

    readReturn = fscanf(file, "%s", uselessString);
    if(readReturn < 1)
    {
        printf("Error reading input file.\n");
        utils_abort();
    }
    readReturn = fscanf(file, "%le", x);
    if(readReturn < 1)
    {
        printf("Error reading input file.\n");
        utils_abort();
    }
}
void utils_readLine(FILE *file, int *x)
{
    int readReturn;
    char uselessString[1024];

    readReturn = fscanf(file, "%s", uselessString);
    if(readReturn < 1)
    {
        printf("Error reading input file.\n");
        utils_abort();
    }
    readReturn = fscanf(file, "%d", x);
    if(readReturn < 1)
    {
        printf("Error reading input file.\n");
        utils_abort();
    }
}
void utils_readLine(FILE *file, char *x)
{
    int readReturn;
    char uselessString[1024];

    readReturn = fscanf(file, "%s", uselessString);
    if(readReturn < 1)
    {
        printf("Error reading input file.\n");
        utils_abort();
    }
    readReturn = fscanf(file, "%s", x);
    if(readReturn < 1)
    {
        printf("Error reading input file.\n");
        utils_abort();
    }
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


