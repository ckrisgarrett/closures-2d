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

#ifdef USE_PAPI
void papi_hwinfo() {
    const PAPI_hw_info_t *hwinfo = NULL;

    if((hwinfo = PAPI_get_hardware_info()) == NULL) {
        printf("PAPI couldn't get hardware info\n");
        utils_abort();
    } else {
        printf("CPU Vendor: %s (%d)\n", hwinfo->vendor_string, hwinfo->vendor);
        printf("CPU Model: %s (%d)\n", hwinfo->model_string, hwinfo->model);
        printf("CPU Revision: %f\n", hwinfo->revision);
        printf("CPUID Family: %d\n", hwinfo->cpuid_family);
        printf("CPUID Model: %d\n", hwinfo->cpuid_model);
        printf("CPUID Stepping: %d\n", hwinfo->cpuid_stepping);
        printf("CPU Max MHz: %d\n", hwinfo->cpu_max_mhz);
        printf("CPU Min MHz: %d\n", hwinfo->cpu_min_mhz);
        printf("Total CPUs: %d\n", hwinfo->totalcpus);
        printf("Sockets: %d\n", hwinfo->sockets);
        printf("Cores per socket: %d\n", hwinfo->cores);
        printf("Hardware threads per core: %d\n", hwinfo->threads);
        printf("Total NUMA Nodes: %d\n", hwinfo->nnodes);
        printf("CPUs per NUMA Node: %d\n", hwinfo->ncpu);
        printf("Virtualized: ");
        if(hwinfo->virtualized)
            printf("yes\n");
        else
            printf("no\n");
        printf("Virtual Vendor: %s\n", hwinfo->virtual_vendor_string);
        printf("Virtual Version: %s\n", hwinfo->virtual_vendor_version);
    }
}

void papi_show_results(int count, int *events, long long *values) {
    PAPI_event_info_t event_info;

    for(int i = 0; i < count; i++) {
        if(PAPI_get_event_info(events[i], &event_info) != PAPI_OK) {
            printf("PAPI invalid event\n");
            utils_abort();
        }
        printf("%s ", event_info.short_descr);
        printf("%lld\n", values[i]);
    }
}
#endif
