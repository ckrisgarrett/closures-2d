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

        sendNorth = (char*)malloc(boundarySizeX);
        recvNorth = (char*)malloc(boundarySizeX);
        sendSouth = (char*)malloc(boundarySizeX);
        recvSouth = (char*)malloc(boundarySizeX);
        sendEast  = (char*)malloc(boundarySizeY);
        recvEast  = (char*)malloc(boundarySizeY);
        sendWest  = (char*)malloc(boundarySizeY);
        recvWest  = (char*)malloc(boundarySizeY);

        if(sendNorth == NULL || sendSouth == NULL || sendEast == NULL || sendWest == NULL ||
           recvNorth == NULL || recvSouth == NULL || recvEast == NULL || recvWest == NULL) {
            printf("comm.cpp: Memory allocation failed.\n");
            utils_abort();
        }
    }
    

    // Get Boundaries.
    getInnerBoundaries(sendNorth, sendSouth, sendEast, sendWest);
    
    
    // Send
        MPI_Isend(sendNorth, boundarySizeX, MPI_CHAR, mpiIndex('y', c_nodeY + 1),
                     sendNorthTag, MPI_COMM_WORLD, &mpiRequest[0]);
        MPI_Isend(sendSouth, boundarySizeX, MPI_CHAR, mpiIndex('y', c_nodeY - 1),
                     sendSouthTag, MPI_COMM_WORLD, &mpiRequest[1]);
        MPI_Isend(sendEast, boundarySizeY, MPI_CHAR, mpiIndex('x', c_nodeX + 1), 
                     sendEastTag, MPI_COMM_WORLD, &mpiRequest[2]);
        MPI_Isend(sendWest, boundarySizeY, MPI_CHAR, mpiIndex('x', c_nodeX - 1), 
                     sendWestTag, MPI_COMM_WORLD, &mpiRequest[3]);
    
    // Receive
    MPI_Irecv(recvNorth, boundarySizeX, MPI_CHAR, mpiIndex('y', c_nodeY + 1),
                     recvNorthTag, MPI_COMM_WORLD, &mpiRequest[4]);
    MPI_Irecv(recvSouth, boundarySizeX, MPI_CHAR, mpiIndex('y', c_nodeY - 1),
                     recvSouthTag, MPI_COMM_WORLD, &mpiRequest[5]);
    MPI_Irecv(recvEast, boundarySizeY, MPI_CHAR, mpiIndex('x', c_nodeX + 1),
                     recvEastTag,MPI_COMM_WORLD, &mpiRequest[6]);
    MPI_Irecv(recvWest, boundarySizeY, MPI_CHAR, mpiIndex('x', c_node - 1),
                     recvWestTag, MPI_COMM_WORLD, &mpiRequest[7]);
    
    // Wait for Send/Recv
    mpiError = MPI_Waitall(8, mpiRequest, MPI_STATUSES_IGNORE);
    if(mpiError != MPI_SUCCESS) {
        printf("comm.cpp: MPI Error %d.\n", mpiError);
        utils_abort();
    }
    
    // Set boundaries.
    setOuterBoundaries(recvNorth, recvSouth, recvEast, recvWest);
}
#else
void Solver::communicateBoundaries()
{
    if(c_initCond == INITCOND_PERIODIC) {
        duplicateBoundaries();
    }
}
#endif
