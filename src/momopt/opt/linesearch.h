/*
 File:   linesearch.h
 Author: Kris Garrett
 Date:   February 13, 2013
*/

#ifndef __LINESEARCH_H
#define __LINESEARCH_H

double linesearch(int n1, int n2, int nq, double *alpha, double *d, double f, double *g, 
                  double *u, double *w, double *p, double ppnDelta, int momentType);

#endif
