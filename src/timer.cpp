#include "timer.h"


Timer::Timer()
{
    timeStarted = false;
    timeEnded = false;
}

void Timer::start()
{
    clock_gettime(CLOCK_MONOTONIC, &time1);
    timeStarted = true;
    timeEnded = false;
}

void Timer::stop()
{
    clock_gettime(CLOCK_MONOTONIC, &time2);
    timeEnded = true;
}

double Timer::getTimeElapsed()
{
    if(!timeStarted || !timeEnded)
        return -1.0;
    
    timespec diffTime;
    diffTime.tv_sec = time2.tv_sec - time1.tv_sec;
    diffTime.tv_nsec = time2.tv_nsec - time1.tv_nsec;
    
    if(diffTime.tv_nsec < 0)
    {
        diffTime.tv_sec--;
        diffTime.tv_nsec += 1E9;
    }
    
    return diffTime.tv_sec + diffTime.tv_nsec / 1E9;
}
