/*
 File:   boundaries.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "kinetic_solver.h"
#include <string.h>


/*
    Get boundary sizes in bytes.
*/
int KineticSolver::getBoundarySizeX()
{
    return c_numQuadPoints * NUM_GHOST_CELLS * (c_gX[3] - c_gX[0] + 1) * sizeof(double);
}

int KineticSolver::getBoundarySizeY()
{
    return c_numQuadPoints * NUM_GHOST_CELLS * (c_gY[3] - c_gY[0] + 1) * sizeof(double);
}


/*
    Returns the inner boundaries for this node's data.
*/
void KineticSolver::getInnerBoundaries(char *north, char *south, char *east, char *west)
{
    int x1 = c_gX[1];
    int x2 = c_gX[2] + 1 - NUM_GHOST_CELLS;
    int y1 = c_gY[1];
    int y2 = c_gY[2] + 1 - NUM_GHOST_CELLS;
    
    
    // North
    if(north != NULL)
    {
        for(int i = c_gX[0]; i <= c_gX[3]; i++)
        {
            for(int j = 0; j < NUM_GHOST_CELLS; j++)
            {
                int index = i + j * (c_gX[3] - c_gX[0] + 1);
                memcpy(&north[index * c_numQuadPoints * sizeof(double)], 
                       &c_kinetic[I3D(i,y2+j,0)], 
                       c_numQuadPoints * sizeof(double));
            }
        }
    }
    
    // South
    if(south != NULL)
    {
        for(int i = c_gX[0]; i <= c_gX[3]; i++)
        {
            for(int j = 0; j < NUM_GHOST_CELLS; j++)
            {
                int index = i + j * (c_gX[3] - c_gX[0] + 1);
                memcpy(&south[index * c_numQuadPoints * sizeof(double)], 
                       &c_kinetic[I3D(i,y1+j,0)], 
                       c_numQuadPoints * sizeof(double));
            }
        }
    }

    // East
    if(east != NULL)
    {
        for(int i = c_gY[0]; i <= c_gY[3]; i++)
        {
            for(int j = 0; j < NUM_GHOST_CELLS; j++)
            {
                int index = i + j * (c_gY[3] - c_gY[0] + 1);
                memcpy(&east[index * c_numQuadPoints * sizeof(double)], 
                       &c_kinetic[I3D(x2+j,i,0)], 
                       c_numQuadPoints * sizeof(double));
            }
        }
    }

    // West
    if(west != NULL)
    {
        for(int i = c_gY[0]; i <= c_gY[3]; i++)
        {
            for(int j = 0; j < NUM_GHOST_CELLS; j++)
            {
                int index = i + j * (c_gY[3] - c_gY[0] + 1);
                memcpy(&west[index * c_numQuadPoints * sizeof(double)],
                       &c_kinetic[I3D(x1+j,i,0)], 
                       c_numQuadPoints * sizeof(double));
            }
        }
    }
}


/*
    Sets the outer boundaries for this node's data.
*/
void KineticSolver::setOuterBoundaries(char *north, char *south, char *east, char *west)
{
    int x0 = c_gX[0];
    int x3 = c_gX[2] + 1;
    int y0 = c_gY[0];
    int y3 = c_gY[2] + 1;
    
    
    // North
    if(north != NULL)
    {
        for(int i = c_gX[0]; i <= c_gX[3]; i++)
        {
            for(int j = 0; j < NUM_GHOST_CELLS; j++)
            {
                int index = i + j * (c_gX[3] - c_gX[0] + 1);
                memcpy(&c_kinetic[I3D(i,y3+j,0)], 
                       &north[index * c_numQuadPoints * sizeof(double)], 
                       c_numQuadPoints * sizeof(double));
            }
        }
    }

    // South
    if(south != NULL)
    {
        for(int i = c_gX[0]; i <= c_gX[3]; i++)
        {
            for(int j = 0; j < NUM_GHOST_CELLS; j++)
            {
                int index = i + j * (c_gX[3] - c_gX[0] + 1);
                memcpy(&c_kinetic[I3D(i,y0+j,0)], 
                       &south[index * c_numQuadPoints * sizeof(double)], 
                       c_numQuadPoints * sizeof(double));
            }
        }
    }

    // East
    if(east != NULL)
    {
        for(int i = c_gY[0]; i <= c_gY[3]; i++)
        {
            for(int j = 0; j < NUM_GHOST_CELLS; j++)
            {
                int index = i + j * (c_gY[3] - c_gY[0] + 1);
                memcpy(&c_kinetic[I3D(x3+j,i,0)], 
                       &east[index * c_numQuadPoints * sizeof(double)], 
                       c_numQuadPoints * sizeof(double));
            }
        }
    }

    // West
    if(west != NULL)
    {
        for(int i = c_gY[0]; i <= c_gY[3]; i++)
        {
            for(int j = 0; j < NUM_GHOST_CELLS; j++)
            {
                int index = i + j * (c_gY[3] - c_gY[0] + 1);
                memcpy(&c_kinetic[I3D(x0+j,i,0)], 
                       &west[index * c_numQuadPoints * sizeof(double)], 
                       c_numQuadPoints * sizeof(double));
            }
        }
    }
}

void KineticSolver::duplicateBoundaries()
{
    int x0 = c_gX[0];
    int x1 = c_gX[1];
    int x2 = c_gX[2] + 1 - NUM_GHOST_CELLS;
    int x3 = c_gX[3] + 1 - NUM_GHOST_CELLS;
    int y0 = c_gY[0];
    int y1 = c_gY[1];
    int y2 = c_gY[2] + 1 - NUM_GHOST_CELLS;
    int y3 = c_gY[3] + 1 - NUM_GHOST_CELLS;

    // x
    for(int i = 0; i < NUM_GHOST_CELLS; i++) {
    for(int j = c_gY[1]; j <= c_gY[2]; j++) {
    for(int q = 0; q < c_numQuadPoints; q++) {
        c_kinetic[I3D(x0+i,j,q)] = c_kinetic[I3D(x2+i,j,q)];
        c_kinetic[I3D(x3+i,j,q)] = c_kinetic[I3D(x1+i,j,q)];
    }}}

    // y
    for(int i = c_gX[1]; i <= c_gX[2]; i++) {
    for(int j = 0; j < NUM_GHOST_CELLS; j++) {
    for(int q = 0; q < c_numQuadPoints; q++) {
        c_kinetic[I3D(i,y0+j,q)] = c_kinetic[I3D(i,y2+j,q)];
        c_kinetic[I3D(i,y3+j,q)] = c_kinetic[I3D(i,y1+j,q)];
    }}}
}
