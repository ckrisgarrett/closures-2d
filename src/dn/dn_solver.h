/*
 File:   dn_solver.h
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#ifndef _DN_SOLVER_H
#define _DN_SOLVER_H

#include "../solver.h"

class DnSolver : public Solver
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
    static const int NUM_GHOST_CELLS = 2;

    void solveFlux(double *moments, double *flux, double dx, double dy);

    int c_numMoments;
    int c_numQuadPoints;
    double *c_moments;
    double *c_flux;
    double *c_w;
    double *c_p;
    double *c_xi;
    double *c_eta;
    
    double *c_fluxOperatorXi;
    double *c_fluxOperatorEta;
    double *c_fluxOperatorXiPlus;
    double *c_fluxOperatorXiMinus;
    double *c_fluxOperatorEtaPlus;
    double *c_fluxOperatorEtaMinus;
    double *c_fluxOperatorXiXi;
    double *c_fluxOperatorXiEta;
    double *c_fluxOperatorEtaEta;
    double *c_flux2OperatorXiXi;
    double *c_flux2OperatorXiEta;
    double *c_flux2OperatorEtaXi;
    double *c_flux2OperatorEtaEta;
};

#endif
