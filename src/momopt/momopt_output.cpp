/*
 File:   output.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "momopt_solver.h"
#include "opt/opt.h"
#include <math.h>
#include <stdio.h>
#include <stdint.h>


/*
    Writes data to filename "out_<time>_<node>.pn".
    Writes optimization stats to "out_<time>_<node>.opt".
    
    For each grid point, the vector of moments is written to the file.  
    The file consists of one long vector.
*/
void MomOptSolver::outputData(double time)
{
    const int FILENAME_SIZE = 256;
    char filename[FILENAME_SIZE];
    FILE *file;
    int64_t temp;
    
    
    // Create and open file.
    snprintf(filename, FILENAME_SIZE, "out_%.3f_%d.pn", time, c_node);
    file = fopen(filename, "wb");
    
    // sizeX and sizeY of output grid
    temp = c_gX[2] - c_gX[1] + 1;
    fwrite(&temp, sizeof(int64_t), 1, file);
    temp = c_gY[2] - c_gY[1] + 1;
    fwrite(&temp, sizeof(int64_t), 1, file);
    
    // Domain bounds
    fwrite(&c_ax, sizeof(double), 1, file);
    fwrite(&c_ay, sizeof(double), 1, file);
    fwrite(&c_bx, sizeof(double), 1, file);
    fwrite(&c_by, sizeof(double), 1, file);
    
    // number of moments
    temp = c_numMoments;
    fwrite(&temp, sizeof(int64_t), 1, file);
    
    // data
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            fwrite(&c_moments[I3D(i,j,0)], sizeof(double), c_numMoments, file);
        }
    }
    
    // Close file.
    fclose(file);
    
    
    // Write optimization data.
    snprintf(filename, FILENAME_SIZE, "out_%.3f_%d.opt", time, c_node);
    
    file = fopen(filename, "w");
    fprintf(file, "Max number of iterations: %d\n", c_optStats.maxIter);
    fprintf(file, "Mean number of iterations: %f\n", c_optStats.iterMean);
    fprintf(file, "Max number of iterations for gamma: %d\n", c_optStats.maxGammaIter);
    fprintf(file, "Mean number of iterations for gamma: %f\n", c_optStats.iterGammaMean);
    fprintf(file, "Number of optimizations: %d\n", c_optStats.numDualSolves);
    fprintf(file, "Regularizations Histogram\n");
    for(int i = 0; i < NUM_REGULARIZATIONS; i++)
    {
        fprintf(file, "r = %.3e   %d\n", REGULARIZATIONS[i], c_optStats.histReg[i]);
    }
    fclose(file);
}

