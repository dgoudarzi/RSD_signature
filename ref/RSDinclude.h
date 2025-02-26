#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

// Choose machine specific macros here or via the compilation command
// define the following if AES instructions are not available and crypto-algorithms is installed (it will be much slower)
// #define GENERICAES   
// define the following if libmd-dev is not installed
// #define GENERICSHA


#ifndef GENERICAES
#include <x86intrin.h> 
#include <wmmintrin.h>  //for intrinsics for AES-NI
#include <tmmintrin.h>

#define aessubkeytype __m128i
#define BENCHMARK_CYCLES
#else
#define aessubkeytype __uint128_t
#include "../crypto-algorithms/aes.h"
#endif

#ifndef GENERICSHA
#include <sha256.h>
#else
#include "../crypto-algorithms/sha256.h"

#define SHA256Init sha256_init
#define SHA256Update sha256_update
#define SHA256Final(hash, ctx_ptr) sha256_final(ctx_ptr, hash)
#endif


#define Lambda 128
#define u64 uint64_t


/* Code Parameters */  
#define NumWindows 217
#define WindowLogModulus 3
#define WindowModulus (1<<WindowLogModulus)
#define CodeLength (WindowModulus*NumWindows)
#define CodeNumChecks 960

#define CodeRepinU128 ((CodeNumChecks+127)>>7)
#define SecretRepinU128 ((4*NumWindows+127)>>7)
#define SecretInnerMask (0x7777777777777777UL)
#define CarryAddMask    (0x8888888888888888UL)
#define SecretLastMask (0x777777777UL)
#define CompressMask (0x7F7F7F7F7F7F7F7FUL)
#define CompressHigh (0x8080808080808080UL)
#define CompressHalf (0x0F0F0F0F0F0F0F0FUL)
#define CompressQuat (0x0303030303030303UL)
#define CompressLow (CompressHigh>>7)
  
#define EXPAND_RATIO 28
#define RSD_EVALUATE

#define AssignNewU64Ptr(TmpPtr, PtrSize) TmpPtr;TmpPtr+=(PtrSize)

#ifdef L8
#define LogRounds 4
#define NbRounds (1L<<LogRounds)
#define LogPlayers 8
#endif

#ifdef L9
#define NbRounds 15
#define LogPlayers 9
#endif

#ifdef L10
#define NbRounds 13
#define LogPlayers 10
#endif

#ifdef L11
#define NbRounds 12
#define LogPlayers 11
#endif

#ifdef L12
#define NbRounds 11
#define LogPlayers 12
#endif

#ifdef L13
#define NbRounds 10
#define LogPlayers 13
#endif

#ifdef L15
#define NbRounds 9
#define LogPlayers 15
#endif

#ifdef L16
#define LogRounds 3
#define NbRounds (1L<<LogRounds)
#define LogPlayers 16
#endif

#define PERM_LEAVES (8*256)

#define NbPlayers (1L<<LogPlayers)
#define ExpRounds (1L<<NbRounds)
#define LEAVES_SIZE (NbPlayers*NbRounds)

#define NB_INNER_RANDOM NbRounds
#define NB_AES_KEYS (2*(LogPlayers+2))

#define uchar unsigned char
#define uint unsigned int

extern void aes128_expand_key(aessubkeytype key_schedule[11], u64 enc_key[2]);
extern void aes128_enc(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule[11]);
extern void aes128_dbl_enc(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule0[11], aessubkeytype key_schedule1[11]);
extern void aes128_dbl_enc_delta_preserve(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule0[11], aessubkeytype key_schedule1[11]);
extern void aes128_expand(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, aessubkeytype key_schedule[11]);
extern void aes128_expand_1cpy(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, aessubkeytype key_schedule[11]);


extern void sha256_dbl_enc(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule0[11], aessubkeytype key_schedule1[11]);
extern void sha256_dbl_enc_delta_preserve(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule0[11], aessubkeytype key_schedule1[11]);
extern void sha256_expand(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, aessubkeytype key_schedule[11]);
extern void sha256_expand_1cpy(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, aessubkeytype key_schedule[11]);


#ifndef DOXYGEN_SHOULD_SKIP_THIS
typedef struct btimer_t {
    unsigned int counter;
    // gettimeofday
    double nb_milliseconds;
    struct timeval start, stop;
    // rdtscp
    uint64_t nb_cycles;
    unsigned int garbage;
    uint64_t cstart, cstop;
} btimer_t;
#endif

extern void btimer_init(btimer_t* timer);
extern void btimer_start(btimer_t *timer);
extern void btimer_count(btimer_t *timer);
extern void btimer_end(btimer_t *timer);
extern double btimer_diff(btimer_t *timer);
extern uint64_t btimer_diff_cycles(btimer_t *timer);
extern double btimer_get(btimer_t *timer);
extern double btimer_get_cycles(btimer_t *timer);
