#include "solver.h"
#include <string.h>
#include <stdlib.h>


/*
 If using MPI, this function communicates the boundary data.
*/
#ifdef USE_MPI
#include <mpi.h>

/*
    Given the node for the x,y direction, returns the 1D index for the node.
*/
int Solver::mpiIndex(char direction, int node)
{
    int nodeX = c_nodeX;
    int nodeY = c_nodeY;
    
    if(direction == 'x')
        nodeX = (node + c_numNodesX) % c_numNodesX;
    else if(direction == 'y')
        nodeY = (node + c_numNodesY) % c_numNodesY;
    else {
        printf("Direction MPI Error.\n");
        utils_abort();
    }
    
    return nodeY + c_numNodesY * nodeX;
}

void Solver::communicateBoundaries()
{
    #ifdef USE_PAPI
    papi_start_update(&c_comm_info);
    #endif
    // Variables.
    int boundarySizeX = getBoundarySizeX();
    int boundarySizeY = getBoundarySizeY();
    int sendNorthTag = 1;
    int recvSouthTag = 1;
    int sendSouthTag = 2;
    int recvNorthTag = 2;
    int sendEastTag  = 3;
    int recvWestTag  = 3;
    int sendWestTag  = 4;
    int recvEastTag  = 4;
    int mpiError;
    MPI_Request mpiRequest[8];
    int requestCount = 0;

    // Variables to allocate only once.
    static bool firstTime = true;
    static char *sendNorth = NULL;
    static char *sendSouth = NULL;
    static char *recvNorth = NULL;
    static char *recvSouth = NULL;
    static char *sendEast = NULL;
    static char *sendWest = NULL;
    static char *recvEast = NULL;
    static char *recvWest = NULL;
    
    if(firstTime)
    {
        firstTime = false;

        if(c_initCond == INITCOND_PERIODIC || c_nodeY < c_numNodesY - 1) {
            sendNorth = (char*)malloc(boundarySizeX);
            recvNorth = (char*)malloc(boundarySizeX);
            if(sendNorth == NULL || recvNorth == NULL) {
                printf("comm.cpp: Memory allocation failed.\n");
                utils_abort();
            }
        }
        if(c_initCond == INITCOND_PERIODIC || c_nodeY > 0) {
            sendSouth = (char*)malloc(boundarySizeX);
            recvSouth = (char*)malloc(boundarySizeX);
            if(sendSouth == NULL || recvSouth == NULL) {
                printf("comm.cpp: Memory allocation failed.\n");
                utils_abort();
            }
        }
        if(c_initCond == INITCOND_PERIODIC || c_nodeX < c_numNodesX - 1) {
            sendEast  = (char*)malloc(boundarySizeY);
            recvEast  = (char*)malloc(boundarySizeY);
            if(sendEast == NULL || recvEast == NULL) {
                printf("comm.cpp: Memory allocation failed.\n");
                utils_abort();
            }
        }
        if(c_initCond == INITCOND_PERIODIC || c_nodeX > 0) {
            sendWest  = (char*)malloc(boundarySizeY);
            recvWest  = (char*)malloc(boundarySizeY);
            if(sendWest == NULL || recvWest == NULL) {
                printf("comm.cpp: Memory allocation failed.\n");
                utils_abort();
            }
        }
    }
    

    // Get Boundaries.
    getInnerBoundaries(sendNorth, sendSouth, sendEast, sendWest);
    
    
    // Send/recv data north.
    if(c_initCond == INITCOND_PERIODIC || c_nodeY < c_numNodesY - 1) {
        MPI_Isend(sendNorth, boundarySizeX, MPI_CHAR, mpiIndex('y', c_nodeY + 1),
                     sendNorthTag, MPI_COMM_WORLD, &mpiRequest[requestCount]);
        MPI_Irecv(recvNorth, boundarySizeX, MPI_CHAR, mpiIndex('y', c_nodeY + 1),
                         recvNorthTag, MPI_COMM_WORLD, &mpiRequest[requestCount + 1]);
        requestCount += 2;
    }

    // Send/recv data south.
    if(c_initCond == INITCOND_PERIODIC || c_nodeY > 0) {
        MPI_Isend(sendSouth, boundarySizeX, MPI_CHAR, mpiIndex('y', c_nodeY - 1),
                     sendSouthTag, MPI_COMM_WORLD, &mpiRequest[requestCount]);
        MPI_Irecv(recvSouth, boundarySizeX, MPI_CHAR, mpiIndex('y', c_nodeY - 1),
                         recvSouthTag, MPI_COMM_WORLD, &mpiRequest[requestCount + 1]);
        requestCount += 2;
    }

    // Send/recv data east.
    if(c_initCond == INITCOND_PERIODIC || c_nodeX < c_numNodesX - 1) {
        MPI_Isend(sendEast, boundarySizeY, MPI_CHAR, mpiIndex('x', c_nodeX + 1), 
                     sendEastTag, MPI_COMM_WORLD, &mpiRequest[requestCount]);
        MPI_Irecv(recvEast, boundarySizeY, MPI_CHAR, mpiIndex('x', c_nodeX + 1),
                         recvEastTag,MPI_COMM_WORLD, &mpiRequest[requestCount + 1]);
        requestCount += 2;
    }

    // Send/recv data west.
    if(c_initCond == INITCOND_PERIODIC || c_nodeX > 0) {
        MPI_Isend(sendWest, boundarySizeY, MPI_CHAR, mpiIndex('x', c_nodeX - 1), 
                     sendWestTag, MPI_COMM_WORLD, &mpiRequest[requestCount]);
        MPI_Irecv(recvWest, boundarySizeY, MPI_CHAR, mpiIndex('x', c_nodeX - 1),
                         recvWestTag, MPI_COMM_WORLD, &mpiRequest[requestCount + 1]);
        requestCount += 2;
    }
    
    // Wait for Send/Recv
    mpiError = MPI_Waitall(requestCount, mpiRequest, MPI_STATUSES_IGNORE);
    if(mpiError != MPI_SUCCESS) {
        printf("comm.cpp: MPI Error %d.\n", mpiError);
        utils_abort();
    }
    
    // Set boundaries.
    setOuterBoundaries(recvNorth, recvSouth, recvEast, recvWest);
    #ifdef USE_PAPI
    papi_finish_update(&c_comm_info);
    #endif
}
#else
void Solver::communicateBoundaries()
{
    #ifdef USE_PAPI
    papi_start_update(&c_comm_info);
    #endif
    if(c_initCond == INITCOND_PERIODIC) {
        duplicateBoundaries();
    }
    #ifdef USE_PAPI
    papi_finish_update(&c_comm_info);
    #endif
}
#endif
