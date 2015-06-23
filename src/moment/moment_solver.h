/*
 File:   moment_solver.h
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#ifndef _MOMENT_SOLVER_H
#define _MOMENT_SOLVER_H

#include "../solver.h"

class MomentSolver : public Solver
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
    enum
    {
        FILTER_TYPE_NONE, FILTER_TYPE_HAUCK, FILTER_TYPE_SSPLINE, FILTER_TYPE_LANCZOS
    };

    static const int NUM_GHOST_CELLS = 2;

    void solveFlux(double *moments, double *flux, double dx, double dy);
    void scaleMoments(double *moments);
    void periodicBoundary(double *u);

    int c_numMoments;
    int c_numQuadPoints;
    double *c_moments;
    double *c_flux;

    double *c_momentFilter;
    int c_filterType;
    double c_filterTune;

    double *c_pnFluxOperatorXiPlus;
    double *c_pnFluxOperatorXiMinus;
    double *c_pnFluxOperatorEtaPlus;
    double *c_pnFluxOperatorEtaMinus;
};

#endif
