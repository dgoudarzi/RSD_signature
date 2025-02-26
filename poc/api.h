#ifndef api_h
#define api_h

#include "pprf_sig.h"

#define CRYPTO_ALGNAME "PPRF SIG"
#define CRYPTO_VERSION "1.00"

#define CRYPTO_SECRETKEYBYTES (((lambda + log_bs*w)/8)+1)
#define CRYPTO_PUBLICKEYBYTES ((small_k/8)+1)
#define CRYPTO_BYTES sizeof(sig_t)

int
crypto_sign_keypair(u64 sk_x[word_x], 
                    u64 seed[word_lambda], 
                    u64 K_0[word_lambda], 
                    u64 K_1[word_lambda], 
                    u64 H[word_H], 
                    u64 y[word_k], 
                    __m128i key_schedule_k0[11], 
                    __m128i key_schedule_k1[11]);

int
crypto_sign(u64 *m, 
            size_t m_len, 
            u64 sk_x[word_x], 
            u64 K_0[word_lambda], 
            u64 K_1[word_lambda], 
            u64 H[word_H], 
            u64 y[word_k], 
            __m128i key_schedule_k0[11], 
            __m128i key_schedule_k1[11],
            sig_t *sigma,
            btimer_t timers_algos[5]
            );

int folding_buckets (__m128i key_schedule_k0[11]);

int
crypto_sign_open(u64 *m, 
                 size_t m_len,
                 u64 H[word_H], 
                 u64 y[word_k], 
                 __m128i key_schedule_k0[11], 
                 __m128i key_schedule_k1[11],
                 sig_t *signature );

#endif /* api_h */