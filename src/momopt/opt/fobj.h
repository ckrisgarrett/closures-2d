/*
 File:   fobj.h
 Author: Kris Garrett
 Date:   May 25, 2012
*/

#ifndef __FOBJ_H
#define __FOBJ_H


// Global variable for Gaunt coefficients.
struct HStructData
{
    int moment;
    double value;
};
struct HStruct
{
    int numMoments;
    HStructData *s;
};
extern HStruct *g_hStruct;


void initFobj(int numOmpThreads, int numManyMoments);
double fobjMn(int n1, int n2, int nq, double *alpha, double *u, double *w, double *p, 
              double *g, double *h, bool useGaunt);
double fobjPPn(int n1, int n2, int nq, double *alpha, double *u, double *w, double *p, double delta, 
               double *g, double *h, bool useGaunt);

#endif
