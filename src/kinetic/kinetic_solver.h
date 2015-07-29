/*
 File:   kinetic_solver.h
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#ifndef _KINETIC_SOLVER_H
#define _KINETIC_SOLVER_H

#include "../solver.h"

class KineticSolver : public Solver
{
public:
    double init(double dx, double dy);
    void update(double dt, double dx, double dy);
    void outputData(double time);
    
    int getBoundarySizeX();
    int getBoundarySizeY();
    void getInnerBoundaries(char *north, char *south, char *east, char *west);
    void setOuterBoundaries(char *north, char *south, char *east, char *west);
    void duplicateBoundaries();

    int getNumGhostCells() { return 2; }

private:
    double *c_kinetic;
    double *c_flux;

    static const int NUM_GHOST_CELLS = 2;

    void initUpdate(int numQuadPoints);
    void solveFlux(double *kinetic, double *flux, double dx, double dy);

    int c_numQuadPoints;

    double *c_xi;
    double *c_eta;
    double *c_quadWeights;
};

#endif
