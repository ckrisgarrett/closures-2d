#include <vector>
#include <string>
#include <map>
#include <string>
#include <papi.h>
#include "profiling.h"
#include "utils.h"

static std::vector<int> profile_events;
static std::map<std::string, PAPI_info_t> profile_info;
static int profile_event_set = PAPI_NULL;

void profile_init() {
    if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
        printf("Couldn't initialize PAPI\n");
        utils_abort();
    }

    if(PAPI_create_eventset(&profile_event_set) != PAPI_OK) {
        printf("PAPI couldn't create event set\n");
        printf("%p\n", &profile_event_set);
        printf("%d\n", PAPI_create_eventset(&profile_event_set));
        utils_abort();
    }

    for(unsigned int i = 0; i < sizeof(default_events) / sizeof(int); i++) {
        if(PAPI_query_event(default_events[i]) == PAPI_OK) {
            if(PAPI_add_event(profile_event_set, default_events[i]) != PAPI_OK) {
                printf("PAPI couldn't add event %x\n", default_events[i]);
            } else {
                profile_events.push_back(default_events[i]);
            }
        }
    }

    if(PAPI_start(profile_event_set) != PAPI_OK) {
        printf("PAPI couldn't start counters\n");
        utils_abort();
    }
}

void profile_start_update(std::string name) {
    profile_info[name]._prev_nsecs = PAPI_get_real_nsec();
    profile_info[name]._prev_cycles = PAPI_get_real_cyc();
    PAPI_read(profile_event_set, profile_info[name]._prev_values);
}

void profile_finish_update(std::string name) {
    long long vals[sizeof(default_events) / sizeof(int)];
    PAPI_read(profile_event_set, vals);

    profile_info[name].cycles += PAPI_get_real_cyc() - profile_info[name]._prev_cycles;
    profile_info[name].nsecs += PAPI_get_real_nsec() - profile_info[name]._prev_nsecs;
    for(unsigned int i = 0; i < profile_events.size(); i++) {
        profile_info[name].values[i] += vals[i] - profile_info[name]._prev_values[i];
    }
    profile_info[name].iterations++;
}

void profile_output() {
    PAPI_event_info_t event_info;

    for(std::map<std::string, PAPI_info_t>::iterator it = profile_info.begin();
            it != profile_info.end(); it++) {
        printf("%s\n", it->first.c_str());
        printf("- Iterations: %lld\n", (it->second).iterations);
        printf("- Nanoseconds: %lld\n", (it->second).nsecs);
        printf("- Cycles: %lld\n", (it->second).cycles);
        for(unsigned int i = 0; i < profile_events.size(); i++) {
            if(PAPI_get_event_info(profile_events[i], &event_info) != PAPI_OK) {
                printf("PAPI invalid event\n");
                utils_abort();
            }
            printf("- %s: ", event_info.short_descr);
            printf("%lld\n", (it->second).values[i]);
        }
    }
}

void profile_stop() {
    if(PAPI_stop(profile_event_set, NULL) != PAPI_OK) {
        printf("PAPI couldn't stop counters\n");
        utils_abort();
    }
}

void profile_hwinfo() {
    const PAPI_hw_info_t *hwinfo = NULL;

    if((hwinfo = PAPI_get_hardware_info()) == NULL) {
        printf("PAPI couldn't get hardware info\n");
        utils_abort();
    } else {
        printf("--- Hardware Info ---\n");
        printf("CPU Vendor: %s (%d)\n", hwinfo->vendor_string, hwinfo->vendor);
        printf("CPU Model: %s (%d)\n", hwinfo->model_string, hwinfo->model);
        printf("CPU Revision: %f\n", hwinfo->revision);
        printf("CPUID Family: %d\n", hwinfo->cpuid_family);
        printf("CPUID Model: %d\n", hwinfo->cpuid_model);
        printf("CPUID Stepping: %d\n", hwinfo->cpuid_stepping);
        printf("CPU Max MHz: %d\n", hwinfo->cpu_max_mhz);
        printf("CPU Min MHz: %d\n", hwinfo->cpu_min_mhz);
        printf("Total CPUs: %d\n", hwinfo->totalcpus);
        printf("Sockets: %d\n", hwinfo->sockets);
        printf("Cores per socket: %d\n", hwinfo->cores);
        printf("Hardware threads per core: %d\n", hwinfo->threads);
        printf("Total NUMA Nodes: %d\n", hwinfo->nnodes);
        printf("CPUs per NUMA Node: %d\n", hwinfo->ncpu);
        printf("Virtualized: ");
        if(hwinfo->virtualized)
            printf("yes\n");
        else
            printf("no\n");
        printf("Virtual Vendor: %s\n", hwinfo->virtual_vendor_string);
        printf("Virtual Version: %s\n", hwinfo->virtual_vendor_version);
    }
}
