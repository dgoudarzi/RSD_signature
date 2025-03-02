#ifndef SDITH_BENCH_UTILS_H
#define SDITH_BENCH_UTILS_H

#include <stdio.h>
#include "api.h"

static inline int get_number_of_tests(int argc, char *argv[], int default_value) {
    int nb_tests = default_value;
    if(argc == 2) {
        if( sscanf(argv[1], "%d", &nb_tests) != 1) {
            printf("Integer awaited.\n");
            return -1;
        } else if( nb_tests <= 0 ) {
            printf("Need to positive integer.\n");
            return -1;
        }
    }
    return nb_tests;
}

static inline void print_configuration() {
    printf("===== SCHEME CONFIG =====\n");
    printf("[API] Algo Name: %s \n", CRYPTO_ALGNAME);
    printf("[PARAM] D:   %d\n", D);
    printf("[PARAM] TAU: %d\n", tau);
}

#endif /* SDITH_BENCH_UTILS_H */
