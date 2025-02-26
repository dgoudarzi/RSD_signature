#ifndef pprf_sig_h
#define pprf_sig_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <openssl/sha.h>
#include <openssl/evp.h>

#include "benchmark/timing.h"

// #define PPRF_SHA


#define u64 uint64_t

#define BS          16
#define big_K       1744       
#define small_k     727
#define lambda      128
#define w           ((big_K)/(BS))
#define D           15
#define tau         9
#define n_par       (1<<D)

#if BS == 16
    #define log_bs 4
#endif

#define word_H (((small_k*big_K)/64)+1)
#define word_K ((big_K/64)+1)
#define word_k ((small_k/64)+1)
#define word_lambda (lambda)/64
#define word_x (((log_bs*w)/64)+1)


typedef struct salt_t {
    u64 K_0[word_lambda];
    u64 K_1[word_lambda];
} salt_t;

typedef struct aux_t {
    u64 x[tau][word_x];
    u64 u_n[tau][word_K];
} aux_t;

typedef struct sig_t {
    salt_t sig_salt;
    unsigned char h1[SHA256_DIGEST_LENGTH];
    unsigned char h2[SHA256_DIGEST_LENGTH];
    u64 copath[tau][D][word_lambda];
    u64 z[tau][word_x];
    aux_t aux;
    u64 com[tau][word_lambda];
    unsigned char com_n[tau][SHA256_DIGEST_LENGTH];
} sig_t;

#include <wmmintrin.h>  //for intrinsics for AES-NI
#include <x86intrin.h>
//compile using gcc and following arguments: -g;-O0;-Wall;-msse2;-msse;-march=native;-maes

//internal stuff

//macros
#define DO_ENC_BLOCK(m,k) \
	do{\
        m = _mm_xor_si128       (m, k[ 0]); \
        m = _mm_aesenc_si128    (m, k[ 1]); \
        m = _mm_aesenc_si128    (m, k[ 2]); \
        m = _mm_aesenc_si128    (m, k[ 3]); \
        m = _mm_aesenc_si128    (m, k[ 4]); \
        m = _mm_aesenc_si128    (m, k[ 5]); \
        m = _mm_aesenc_si128    (m, k[ 6]); \
        m = _mm_aesenc_si128    (m, k[ 7]); \
        m = _mm_aesenc_si128    (m, k[ 8]); \
        m = _mm_aesenc_si128    (m, k[ 9]); \
        m = _mm_aesenclast_si128(m, k[10]);\
    }while(0)

#define DO_DEC_BLOCK(m,k) \
	do{\
        m = _mm_xor_si128       (m, k[10+0]); \
        m = _mm_aesdec_si128    (m, k[10+1]); \
        m = _mm_aesdec_si128    (m, k[10+2]); \
        m = _mm_aesdec_si128    (m, k[10+3]); \
        m = _mm_aesdec_si128    (m, k[10+4]); \
        m = _mm_aesdec_si128    (m, k[10+5]); \
        m = _mm_aesdec_si128    (m, k[10+6]); \
        m = _mm_aesdec_si128    (m, k[10+7]); \
        m = _mm_aesdec_si128    (m, k[10+8]); \
        m = _mm_aesdec_si128    (m, k[10+9]); \
        m = _mm_aesdeclast_si128(m, k[0]);\
    }while(0)

#define AES_128_key_exp(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))

#endif
