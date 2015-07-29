#ifndef _PROFILE_H
#define _PROFILE_H

#include <string>
#include <papi.h>

const int default_events[] = {
    PAPI_L1_DCH,
    PAPI_L1_DCM,
    PAPI_L1_ICH,
    PAPI_L1_ICM,
    PAPI_L2_DCH,
    PAPI_L2_DCM,
    PAPI_L2_ICH,
    PAPI_L2_ICM,
    PAPI_L3_DCH,
    PAPI_L3_DCM,
    PAPI_L3_ICH,
    PAPI_L3_ICM,
    PAPI_INT_INS,
    PAPI_FP_INS,
    PAPI_FP_OPS,
    PAPI_SP_OPS,
    PAPI_DP_OPS,
    PAPI_VEC_INS,
    PAPI_VEC_SP,
    PAPI_VEC_DP
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
