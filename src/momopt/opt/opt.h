/*
 File:   opt.h
 Author: Kris Garrett
 Date:   February 13, 2013
*/

#ifndef __OPT_H
#define __OPT_H


struct OPTIONS
{
    int momentType;
    int maxIter;
    int maxBfgsIter;        // For BFGS only
    double tolAbs;
    double tolRel;
    double tolGamma;
    double condHMax;
    double condHMaxBfgs;    // For BFGS only
    double deltaPPn;        // For PPn only
    bool useGaunt;  // For Mn non-change of basis only
};

struct OUTS
{
    double r;
    int iter;
    int iterGamma;
    double gamma;
    double normG;
};

enum
{
    MOMENT_TYPE_MN, MOMENT_TYPE_PPN
};

static const int NUM_REGULARIZATIONS = 5;
static const double REGULARIZATIONS[NUM_REGULARIZATIONS] = {0.0, 1e-3, 1e-2, 1e-1, 5e-1};


void opt(int nm, int nm2, int nq, double *u, double *alpha, OPTIONS options, double *w, double *p, 
         OUTS *outs);
void optcb(int nm, int nq, double *u, double *P2M, double *alphaP, double *alphaM, 
           double *w, double *p, OPTIONS options, OUTS *outs);
void optbfgs(int nm, int nq, double *u, double *P2M, double *alphaP, double *alphaM, 
             double *w, double *p, OPTIONS options, OUTS *outs);

#endif
