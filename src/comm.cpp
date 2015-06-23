#include "solver.h"
#include <string.h>
#include <stdlib.h>


/*
 If using MPI, this function communicates the boundary data.
*/
#ifdef USE_MPI
#include <mpi.h>

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
    MPI_Status  mpiStatus;

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
        if(c_nodeY < c_numNodesY - 1)
        {
            sendNorth = (char*)malloc(boundarySizeX);
            recvNorth = (char*)malloc(boundarySizeX);
        }
        if(c_nodeY > 0)
        {
            sendSouth = (char*)malloc(boundarySizeX);
            recvSouth = (char*)malloc(boundarySizeX);
        }
        if(c_nodeX < c_numNodesX - 1)
        {
            sendEast  = (char*)malloc(boundarySizeY);
            recvEast  = (char*)malloc(boundarySizeY);
        }
        if(c_nodeX > 0)
        {
            sendWest  = (char*)malloc(boundarySizeY);
            recvWest  = (char*)malloc(boundarySizeY);
        }
    }
    

    // Get Boundaries.
    getInnerBoundaries(sendNorth, sendSouth, sendEast, sendWest);
    
    
    // Send/Recv data north.
    if(c_nodeY < c_numNodesY - 1)
    {
        MPI_Sendrecv(sendNorth, boundarySizeX, MPI_CHAR, c_node + 1, 
                     sendNorthTag, 
                     recvNorth, boundarySizeX, MPI_CHAR, c_node + 1,
                     recvNorthTag, 
                     MPI_COMM_WORLD, &mpiStatus);
    }
    
    // Send data south.
    if(c_nodeY > 0)
    {
        MPI_Sendrecv(sendSouth, boundarySizeX, MPI_CHAR, c_node - 1, 
                     sendSouthTag, 
                     recvSouth, boundarySizeX, MPI_CHAR, c_node - 1,
                     recvSouthTag,
                     MPI_COMM_WORLD, &mpiStatus);
    }
    
    // Send data east.
    if(c_nodeX < c_numNodesX - 1)
    {
        MPI_Sendrecv(sendEast, boundarySizeY, MPI_CHAR, c_node + c_numNodesY, 
                     sendEastTag, 
                     recvEast, boundarySizeY, MPI_CHAR, c_node + c_numNodesY,
                     recvEastTag,
                     MPI_COMM_WORLD, &mpiStatus);
    }
    
    // Send data west.
    if(c_nodeX > 0)
    {
        MPI_Sendrecv(sendWest, boundarySizeY, MPI_CHAR, c_node - c_numNodesY, 
                     sendWestTag, 
                     recvWest, boundarySizeY, MPI_CHAR, c_node - c_numNodesY,
                     recvWestTag,
                     MPI_COMM_WORLD, &mpiStatus);
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
