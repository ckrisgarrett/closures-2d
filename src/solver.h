/*
 File:   solver.h
 Author: Kris Garrett
 Date:   February, 2013
*/

#ifndef _SOLVER_H
#define _SOLVER_H

#define I2D(i,j) ((i)*(c_gY[3]-c_gY[0]+1) + (j))
#define I3D(i,j,k) (((i)*(c_gY[3]-c_gY[0]+1) + (j))*(c_vectorSize) + (k))

class Solver
{
public:
    virtual double init(double dx, double dy) = 0;
    virtual int getNumGhostCells() = 0;
    virtual void update(double dt, double dx, double dy) = 0;
    virtual void outputData(double time) = 0;
    
    virtual int getBoundarySizeX() = 0;
    virtual int getBoundarySizeY() = 0;
    virtual void getInnerBoundaries(char *north, char *south, char *east, char *west) = 0;
    virtual void setOuterBoundaries(char *north, char *south, char *east, char *west) = 0;
    virtual void duplicateBoundaries() = 0;
    
    void communicateBoundaries();
    void initializeGrid(int numGhostCells, double dx, double dy, double sigma);
    
    int c_node;                     // MPI node index
    int c_numNodes;                 // Number of MPI nodes
    int c_nodeX, c_nodeY;           // Node index for x,y-direction
    int c_numNodesX, c_numNodesY;   // Number of nodes for the x,y-direction
    int c_numOmpThreads;            // Number of openmp threads
    
    double *c_initialGrid;          // Initial grid
    double *c_sigmaT;               // sigma_t at each point of the grid
    double *c_sigmaS;               // sigma_s at each point of the grid
    int c_sizeX, c_sizeY;           // Number of points in the x,y-direction
    double c_ax, c_ay, c_bx, c_by;  // Rectangular domain bounds [ax,ay] x [bx,by]
    double c_initX, c_initY;        // Least x,y value on the grid
    
    int c_gX[4];                    // Indices for ghost values in x-direction
    int c_gY[4];                    // Indices for ghost values in y-direction
    
    double c_gaussianSigma;         // Standard deviation for guassian initial condition
    double c_floor;                 // Smallest value for the density initial condition
    int c_initCond;                 // 0: Gaussian/Linesource   1: Lattice   2: Periodic
    int c_vectorSize;               // number of moments or quadrature points (Sn)
    
    enum
    {
        INITCOND_LINESOURCE, INITCOND_LATTICE, INITCOND_PERIODIC
    };
};

#endif
