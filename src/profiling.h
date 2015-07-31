#ifndef _PROFILE_H
#define _PROFILE_H

#include <string>
#include <papi.h>

const int default_events[] = {
    PAPI_L1_TCH,
    PAPI_L1_TCM,
    PAPI_L2_TCH,
    PAPI_L2_TCM,
    PAPI_L3_TCH,
    PAPI_L3_TCM,
};

typedef struct _PAPI_info {
    long long iterations;
    long long cycles;
    long long nsecs;
    long long values[sizeof(default_events) / sizeof(int)];
    long long _prev_cycles;
    long long _prev_nsecs;
    long long _prev_values[sizeof(default_events) / sizeof(int)];
} PAPI_info_t;

void profile_start_update(std::string);
void profile_finish_update(std::string);
void profile_show_results();
void profile_stop();
void profile_init();
void profile_output();
void profile_hwinfo();
#endif
