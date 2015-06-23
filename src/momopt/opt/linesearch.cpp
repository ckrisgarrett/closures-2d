/*
 File:   linesearch.cpp
 Author: Kris Garrett
 Date:   February 13, 2013
*/

#include "linesearch.h"
#include "fobj.h"
#include "opt.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>


/*
    Find alpha which decreases f.
*/
double linesearch(int n1, int n2, int nq, double *alpha, double *d, double f, double *g, 
                  double *u, double *w, double *p, double ppnDelta, int momentType)
{
    double t = 1.0;
    double fNew = f;
    double *alphaOld = (double*)malloc(n1 * sizeof(double));
    
    for(int i = 0; i < n1; i++)
        alphaOld[i] = alpha[i];
    
    
    // fixed sufficient decrease parameter (this parameter is
    // often denoted alpha in literature discussing Armijo's method
    double delta = 1e-3; 
    double beta = .5;       // decrease rate of step size
    
    
    // Create factor to check for backtracking.
    double maxFactor = 0;
    for(int i = 0; i < n1; i++)
    {
        double temp = fabs(d[i]) / (fabs(alphaOld[i]) + 1e-100);
        if(temp > maxFactor)
            maxFactor = temp;
    }
    
    
    // backtracking
    while(t * maxFactor > DBL_EPSILON / 2)
    {
        // alpha = alphaOld + t*d;
        for(int i = 0; i < n1; i++)
            alpha[i] = alphaOld[i] + t * d[i];
        
        switch(momentType)
        {
            case MOMENT_TYPE_MN:
                fNew = fobjMn(n1, n2, nq, alpha, u, w, p, NULL, NULL, false);
                break;
            case MOMENT_TYPE_PPN:
                fNew = fobjPPn(n1, n2, nq, alpha, u, w, p, ppnDelta, NULL, NULL, false);
                break;
            default:
                break;
        }
        
        // check the Armijo criterion with the current quadrature
        double temp = 0.0;
        for(int i = 0; i < n1; i++)
            temp += g[i] * d[i];
        if(fNew <= f + delta * t * temp)
        {
            free(alphaOld);
            return t;
        }
        
        t = beta * t;
    }
    
    
    // Getting to here means we backtracked all the way.
    // if we backtracked all the way back to the starting point, here we make
    // sure to send back evaluations of the objective function and gradient
    // using the (possibly) new quadrature
    for(int i = 0; i < n1; i++)
        alpha[i] = alphaOld[i];
    free(alphaOld);
    return 0.0;
}


