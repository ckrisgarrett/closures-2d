/*
 File:   opt.cpp
 Author: Kris Garrett
 Date:   February 13, 2013
*/

#include "opt.h"
#include "fobj.h"
#include "linesearch.h"
#include "../../utils.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>


static
void optScaled(int nm, int nm2, int nq, double *u, double *alpha, double *w, double *p, 
               OPTIONS *options, OUTS *outs);


/*
    Calculate alpha for u = (1,0,0,...,0)^T
*/
static void isotropicAlpha(int n, double *alpha, int momentType, double delta)
{
    switch(momentType)
    {
        case MOMENT_TYPE_MN:
            alpha[0] = 2.0 * sqrt(M_PI) * log(0.5 / sqrt(M_PI));
            for(int i = 1; i < n; i++)
                alpha[i] = 0.0;
            break;
        case MOMENT_TYPE_PPN:
            alpha[0] = 1.0 - 4.0 * delta * M_PI;
            for(int i = 1; i < n; i++)
                alpha[i] = 0.0;
            break;
        default:
            printf("isotropicAlpha: moment type out of range.\n");
            utils_abort();
            break;
    }
}


/*
    Estimate an upper bound for the error.
*/
static double modcfl(int n, double *da, int momentType)
{
    double da_1 = 0;
    for(int i = 0; i < n; i++)
    {
        da_1 = da_1 + fabs(da[i]);
    }
    
    int pnOrder = floor(sqrt(n)) + 1;
    double m_inf = sqrt((2 * pnOrder + 1) / (2 * M_PI));
    
    switch(momentType)
    {
        case MOMENT_TYPE_MN:
            return exp(2.0 * da_1 * m_inf);
        case MOMENT_TYPE_PPN:
            return 1.0 + 2.0 * da_1 * m_inf;
        default:
            printf("modcfl: moment type out of range.\n");
            utils_abort();
    }
    
    // Should never reach this point.
    return 0.0;
}


/*
    Scales u and calls the optimization routine.
*/
void opt(int nm, int nm2, int nq, double *u, double *alpha, OPTIONS options, double *w, double *p, 
         OUTS *outs)
{
    switch(options.momentType)
    {
        case MOMENT_TYPE_MN:
        {
            // Scale u and alpha.
            double u0 = u[0];
            double alphaScale = log(u0) * 2.0 * sqrt(M_PI);
            alpha[0] = alpha[0] - alphaScale;
            
            for(int i = 0; i < nm; i++)
                u[i] = u[i] / u0;
            
            if(u0 <= 0)
            {
                printf("u[0] < 0\n");
                utils_abort();
            }
            
            // Run optimization.
            optScaled(nm, nm2, nq, u, alpha, w, p, &options, outs);

            // Rescale u and alpha.
            for(int i = 0; i < nm; i++)
                u[i] = u[i] * u0;
            alpha[0] = alpha[0] + alphaScale;

            break;
        }
        case MOMENT_TYPE_PPN:
        {
            // Scale u and alpha.
            double u0 = u[0];
            for(int i = 0; i < nm; i++)
            {
                alpha[i] = alpha[i] / u0;
                u[i] = u[i] / u0;
            }
            options.deltaPPn = options.deltaPPn / (u0 * u0);
            
            if(u0 <= 0)
            {
                printf("u[0] < 0\n");
                utils_abort();
            }
    
            // Run optimization.
            optScaled(nm, nm2, nq, u, alpha, w, p, &options, outs);
            
            // Rescale u and alpha.
            for(int i = 0; i < nm; i++)
            {
                alpha[i] = alpha[i] * u0;
                u[i] = u[i] * u0;
            }

            break;
        }
        default:
        {
            printf("Moment type not supported.\n");
            utils_abort();
        }
    }
}


/*
    Does the optimization for u = (1,...)^T
*/
void optScaled(int nm, int nm2, int nq, double *u, double *alpha, double *w, double *p, 
               OPTIONS *options, OUTS *outs)
{
    double f = 0;
    double *d = (double*)malloc(nm * sizeof(double));
    double *cholesky = (double*)malloc(nm * nm * sizeof(double));
    double *g = (double*)malloc(nm * sizeof(double));
    double *h = (double*)malloc(nm * nm * sizeof(double));
    double *uOriginal = (double*)malloc(nm * sizeof(double));
    double *alpha0 = (double*)malloc(nm * sizeof(double));
    double u0;
    
    int maxIter = options->maxIter;
    double tolRel = options->tolRel;
    double tolGamma = options->tolGamma;
    double r = 0.0;
    double condHMax = options->condHMax;
    int regIndex = 0;
    int inc1 = 1;
    char lChar = 'L';
    
    for(int i = 0; i < nm; i++)
        uOriginal[i] = u[i];
    
    
    // tolerance relative to |u|
    double rtol = dnrm2_(&nm, u, &inc1) * tolRel;
    
    
    // initialize the iteration counter
    int iter = 0;
    int totalIter = 0;
    
    // initialize a variable to hold the number of steps we took get get gamma
    // below its tolerance after the norm of the gradient satisfies its
    // tolerance.  (typically, the norm of the gradient will fall below its
    // tolerance much more quickly than gamma will.)
    int iterGamma = 0;

    // initialize the backtracking factor
    double t = 1;
    double gamma, err;
    gamma = err = 0;      // To stop compiler warnings.
    
    // set the multipliers we'll use to iterate to the initial condition
    for(int i = 0; i < nm; i++)
        alpha0[i] = alpha[i];
    
    
    // Iterate to find alpha.
    while(true)
    {
        // increment the iteration counters
        iter = iter + 1;
        totalIter = totalIter + 1;

        // Find f, g, and h.
        switch(options->momentType)
        {
            case MOMENT_TYPE_MN:
                f = fobjMn(nm, nm2, nq, alpha, u, w, p, g, h, options->useClebschGordan);
                break;
            case MOMENT_TYPE_PPN:
                f = fobjPPn(nm, nm2, nq, alpha, u, w, p, options->deltaPPn, g, h, 
                    options->useClebschGordan);
        }

        
        // Check quadrature and find Newton direction with Cholesky factorization.
        // ALSO NEED TO CHECK CONDITION NUMBER OF H.
        memcpy(cholesky, h, nm * nm * sizeof(double));
        int error;
        dpotrf_(&lChar, &nm, cholesky, &nm, &error);
        double norm = utils_norm1(nm, h);
        double condest = 0.0;
        int error2;
        double *work = (double*)malloc(3*nm*sizeof(double));
        int *iwork = (int*)malloc(nm*sizeof(int));
        dpocon_(&lChar, &nm, cholesky, &nm, &norm, &condest, work, iwork, &error2);
        free(work);
        free(iwork);
        condest = 1.0 / condest;


        // Initialize to isotropic alpha if H is too bad on first iteration.
//        if(totalIter == 1 && (error != 0 || condest > condHMax))
//        {
//            isotropicAlpha(nm, alpha0, options->momentType, options->deltaPPn);
//            
//            // Reset Variables.
//            rtol = dnrm2_(&nm, u, &inc1) * tolRel;
//            memcpy(alpha, alpha0, nm * sizeof(double));
//            iter = 0;
//            iterGamma = 0;
//            t = 1;
//            
//            continue;
//        }

        
        // Regularize criteria.
        if(iter >= maxIter || error != 0 || t == 0 || condest > condHMax)
        {
            // REGULARIZE ...
            if(regIndex < NUM_REGULARIZATIONS - 1)
            {
                regIndex++;
                r = REGULARIZATIONS[regIndex];
                
                u[0] = 1.0;
                for(int i = 1; i < nm; i++)
                    u[i] = uOriginal[i] * (1.0 - r);
            }
            else
            {
                // This should never happen.
                printf("Max regularization achieved.\n");
                printf("%d: %d: %e: %e\n", iter, error, t, condest);
                utils_abort();
            }
            
            // TEST CODE
            isotropicAlpha(nm, alpha0, options->momentType, options->deltaPPn);
            // TEST CODE
            
            // Reset variables.
            rtol = dnrm2_(&nm, u, &inc1) * tolRel;
            memcpy(alpha, alpha0, nm * sizeof(double));
            iterGamma = 0;
            iter = 0;
            t = 1;
            if(regIndex == NUM_REGULARIZATIONS - 1)
                maxIter = 10 * maxIter;
            
            continue;
        }
        
        
        // Solve for line direction.
        for(int i = 0; i < nm; i++)
            d[i] = g[i];
        dpotrs_(&lChar, &nm, &inc1, cholesky, &nm, d, &nm, &error);
        for(int i = 0; i < nm; i++)
            d[i] = -d[i];
        

        // the error as measured by the norm of the gradient
        err = dnrm2_(&nm, g, &inc1);
        
        // gamma, the supplementary stopping criterion
        gamma = modcfl(nm, d, options->momentType);
        
        
        // now check the stopping criteria
        if(err <= rtol && gamma <= tolGamma)
        {
            // problem solved!
            if(options->momentType == MOMENT_TYPE_MN) {
                u0 = u[0] + g[0];
                for(int i = 0; i < nm; i++) {
                    u[i] += g[i];
                    u[i] /= u0;
                }
                alpha[0] -= log(u0) * 2.0 * sqrt(M_PI);
            }

            break;
        }
        else if(err <= rtol && gamma > tolGamma)
        {
            // this is to keep track of how many extra iterations we took just for gamma
            iterGamma++;
        }
        
        
        // if you get to this point, that means you didn't satisfy the stopping
        // criteria and you didn't need to regularize more than you already are
        // so, it's time to take a step
        t = linesearch(nm, nm2, nq, alpha, d, f, g, u, w, p, options->deltaPPn, options->momentType);
    }
    
    
    if(outs != NULL)
    {
        outs->iter = totalIter;
        outs->iterGamma = iterGamma;
        if(outs->iterGamma < 0)
            outs->iterGamma = 0;
        outs->gamma = gamma;
        outs->normG = err;
        outs->r = r;
    }

    
    free(d);
    free(cholesky);
    free(g);
    free(h);
    free(uOriginal);
    free(alpha0);
}


