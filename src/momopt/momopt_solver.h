/*
 File:   momopt_solver.h
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#ifndef _MOMOPT_SOLVER_H
#define _MOMOPT_SOLVER_H

#include "../solver.h"
#include "opt/opt.h"

#define I3D2(i,j,k) (((i)*(c_gY[3]-c_gY[0]+1) + (j))*(c_numMoments*c_numMoments) + (k))

class MomOptSolver : public Solver
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
    
    struct OPTIMIZATION_STATS
    {
        int maxIter;
        int maxGammaIter;
        int histReg[NUM_REGULARIZATIONS];
        double iterMean;
        double iterGammaMean;
        int numDualSolves;
    };

    enum
    {
        OPT_TYPE_ORIGINAL, OPT_TYPE_CHANGE_BASIS, OPT_TYPE_BFGS
    };

    
    void initUpdate(int numMoments);
    void solveFlux(double *moments, double *flux, double *alpha, double dx, double dy);
    void solveOptimization();
    double computeAnsatz(double *alpha, int q);

    int c_numMoments;
    int c_numManyMoments;
    int c_numQuadPoints;

    int c_optType;
    int c_momentType;
    double c_deltaPPn;
    double c_condHMax;
    double c_condHMaxBfgs;
    int c_useGaunt;
    double c_tolerance;
    int c_maxIter;
    int c_maxIterBfgs;
    double c_theta;
    
    double *c_spHarm;
    double *c_quadWeights;
    double *c_xi;
    double *c_eta;
    
    double *c_moments;
    double *c_alpha;
    double *c_flux;
    double *c_alphaP;
    double *c_P2M;
    
    OPTIMIZATION_STATS c_optStats;
};

#endif
