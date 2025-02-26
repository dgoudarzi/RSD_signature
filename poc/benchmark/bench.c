#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "timing.h"
#include "utils.h"
#include "api.h"

#define B_KEY_GENERATION 0
#define B_SIGN_ALGO 1
#define B_VERIFY_ALGO 2
#define B_TREE_GEN 3
#define B_FOLDING_BUCKETS 4
#define NUMBER_OF_ALGO_BENCHES 5


int main(int argc, char *argv[]) {
    srand((unsigned int) time(NULL));

    int nb_tests = get_number_of_tests(argc, argv, 1);
    if(nb_tests < 0)
        exit(EXIT_FAILURE);

    print_configuration();
    printf("\n");

    btimer_t timers_algos[NUMBER_OF_ALGO_BENCHES];
    double std_timer[NUMBER_OF_ALGO_BENCHES];

    // Initialisation
    double timer_pow2[NUMBER_OF_ALGO_BENCHES];
    for(int j=0; j<NUMBER_OF_ALGO_BENCHES; j++) {
        btimer_init(&timers_algos[j]);
        timer_pow2[j] = 0;
    }

    // Execution
    int score = 0;
    int ret;
    u64 sk_x[word_x];
    u64 sk_seed[word_lambda] = {0};
    u64 K_0[word_lambda], K_1[word_lambda];
    __m128i ks_k0[11], ks_k1[11];
    u64 H[word_H];
    u64 y[word_k];
    sig_t sig;
    for(int i=0; i<nb_tests; i++) {
        // KEYGEN
        for (size_t i=0; i<word_x; i++) {
            sk_x[i] = 0;
        }
        btimer_start(&timers_algos[B_KEY_GENERATION]);
        ret = crypto_sign_keypair(sk_x, sk_seed, K_0, K_1, H, y, ks_k0, ks_k1);
        btimer_end(&timers_algos[B_KEY_GENERATION]);
        btimer_count(&timers_algos[B_KEY_GENERATION]);
        timer_pow2[B_KEY_GENERATION] += pow(btimer_diff(&timers_algos[B_KEY_GENERATION]), 2) / nb_tests;
        if(ret) {
            printf("Failure (num %d): crypto_sign_keypair\n", i);
            continue;
        }

        // Select the message
        #define MLEN 32
        uint64_t m[MLEN] = {1, 2, 3, 4};
        
        // Sign the message
        btimer_start(&timers_algos[B_SIGN_ALGO]);
        ret = crypto_sign(m, MLEN, sk_x, K_0, K_1, H, y, ks_k0, ks_k1, &sig, timers_algos);
        btimer_end(&timers_algos[B_SIGN_ALGO]);
        btimer_count(&timers_algos[B_SIGN_ALGO]);
        timer_pow2[B_SIGN_ALGO] += pow(btimer_diff(&timers_algos[B_SIGN_ALGO]), 2) / nb_tests;
        btimer_count(&timers_algos[B_TREE_GEN]);
        timer_pow2[B_TREE_GEN] += pow(btimer_diff(&timers_algos[B_TREE_GEN]), 2) / nb_tests;

        btimer_start(&timers_algos[B_FOLDING_BUCKETS]);
        ret = folding_buckets(ks_k0);
        btimer_end(&timers_algos[B_FOLDING_BUCKETS]);
        btimer_count(&timers_algos[B_FOLDING_BUCKETS]);
        timer_pow2[B_FOLDING_BUCKETS] += pow(btimer_diff(&timers_algos[B_FOLDING_BUCKETS]), 2) / nb_tests;
        
        if(ret) {
            printf("Failure (num %d): crypto_sign\n", i);
            continue;
        }

        // Verify/Open the signature
        btimer_start(&timers_algos[B_VERIFY_ALGO]);
        ret = crypto_sign_open(m, MLEN, H, y, ks_k0, ks_k1, &sig);
        btimer_end(&timers_algos[B_VERIFY_ALGO]); btimer_count(&timers_algos[B_VERIFY_ALGO]);
        timer_pow2[B_VERIFY_ALGO] += pow(btimer_diff(&timers_algos[B_VERIFY_ALGO]), 2) / nb_tests;
        if(ret) {
            printf("Failure (num %d): crypto_sign_open\n", i);
            continue;
        }
        
        score++;
    }

    // Compute some statistics
    std_timer[B_KEY_GENERATION] = sqrt(timer_pow2[B_KEY_GENERATION] - pow(btimer_get(&timers_algos[B_KEY_GENERATION]),2));
    std_timer[B_SIGN_ALGO] = sqrt(timer_pow2[B_SIGN_ALGO] - pow(btimer_get(&timers_algos[B_SIGN_ALGO]),2));
    std_timer[B_VERIFY_ALGO] = sqrt(timer_pow2[B_VERIFY_ALGO] - pow(btimer_get(&timers_algos[B_VERIFY_ALGO]),2));
    std_timer[B_TREE_GEN] = sqrt(timer_pow2[B_TREE_GEN] - pow(btimer_get(&timers_algos[B_TREE_GEN]),2));
    std_timer[B_FOLDING_BUCKETS] = sqrt(timer_pow2[B_FOLDING_BUCKETS] - pow(btimer_get(&timers_algos[B_FOLDING_BUCKETS]),2));

    // Display Infos
    printf("===== SUMMARY =====\n");
    printf("Correctness: %d/%d\n", score, nb_tests);
    printf("\n");

    printf("Timing in ms:\n");
    printf(" - Key Gen: %.2f ms (std=%.2f)\n",
        btimer_get(&timers_algos[B_KEY_GENERATION]),
        std_timer[B_KEY_GENERATION]
    );
    printf(" - Sign:    %.2f ms (std=%.2f)\n",
        btimer_get(&timers_algos[B_SIGN_ALGO]),
        std_timer[B_SIGN_ALGO]
    );
    printf(" - Verify:  %.2f ms (std=%.2f)\n",
        btimer_get(&timers_algos[B_VERIFY_ALGO]),
        std_timer[B_VERIFY_ALGO]
    );
    printf(" - TreeGen:  %.2f ms (std=%.2f)\n",
        btimer_get(&timers_algos[B_TREE_GEN]),
        std_timer[B_TREE_GEN]
    );
    printf(" - FoldingBuckets:  %.2f ms (std=%.2f)\n",
        btimer_get(&timers_algos[B_FOLDING_BUCKETS]),
        std_timer[B_FOLDING_BUCKETS]
    );
    printf("\n");

    printf("Timing in cycles:\n");
    printf(" - Key Gen: %.2f cycles\n", btimer_get_cycles(&timers_algos[B_KEY_GENERATION]));
    printf(" - Sign:    %.2f cycles\n", btimer_get_cycles(&timers_algos[B_SIGN_ALGO]));
    printf(" - Verify:  %.2f cycles\n", btimer_get_cycles(&timers_algos[B_VERIFY_ALGO]));
    printf(" - TreeGen:  %.2f cycles\n", btimer_get_cycles(&timers_algos[B_TREE_GEN]));
    printf(" - FoldingBucket:  %.2f cycles\n", btimer_get_cycles(&timers_algos[B_FOLDING_BUCKETS]));
    
    printf("\n");


    printf("Communication cost:\n");
    printf(" - PK size: %d B\n", CRYPTO_PUBLICKEYBYTES);
    printf(" - SK size: %d B\n", CRYPTO_SECRETKEYBYTES);
    printf(" - Signature size (MAX): %ld B\n", CRYPTO_BYTES);
    printf("\n");

    return 0;
}
