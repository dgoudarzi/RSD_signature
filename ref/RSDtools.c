/**
 * @file
 * This file contains low-level tools necessary for RSD signatures.
 **/

#include "RSDinclude.h"

#ifndef GENERICAES
#ifndef DOXYGEN_SHOULD_SKIP_THIS
static __m128i aes_128_key_expansion(__m128i key, __m128i keygened){
    keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    return _mm_xor_si128(key, keygened);
}

#define AES_128_key_exp_call_round(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))
#endif

/*! \brief Key expansion of an AES-128 key into 11 round subkeys

Note: All the AES-NI code comes from or is inspired by https://gist.github.com/acapola/d5b940da024080dfaf5f
*/
void aes128_expand_key(__m128i key_schedule[11], u64 enc_key[2]){
    key_schedule[0] = _mm_set_epi64x(enc_key[0], enc_key[1]);
    key_schedule[1]  = AES_128_key_exp_call_round(key_schedule[0], 0x01);
    key_schedule[2]  = AES_128_key_exp_call_round(key_schedule[1], 0x02);
    key_schedule[3]  = AES_128_key_exp_call_round(key_schedule[2], 0x04);
    key_schedule[4]  = AES_128_key_exp_call_round(key_schedule[3], 0x08);
    key_schedule[5]  = AES_128_key_exp_call_round(key_schedule[4], 0x10);
    key_schedule[6]  = AES_128_key_exp_call_round(key_schedule[5], 0x20);
    key_schedule[7]  = AES_128_key_exp_call_round(key_schedule[6], 0x40);
    key_schedule[8]  = AES_128_key_exp_call_round(key_schedule[7], 0x80);
    key_schedule[9]  = AES_128_key_exp_call_round(key_schedule[8], 0x1B);
    key_schedule[10] = AES_128_key_exp_call_round(key_schedule[9], 0x36);
}

/*! \brief Standard AES-ECB mode provided by https://gist.github.com/acapola/d5b940da024080dfaf5f */
void aes128_enc(u64 *plainText, u64 *cipherText, u64 nbblocks, __m128i key_schedule[11]){
  __m128i m;

  for(u64 i=0;i<nbblocks;i++,plainText+=2,cipherText+=2) {
    m = _mm_loadu_si128((__m128i *) plainText);
    m = _mm_xor_si128       (m, key_schedule[ 0]);
    m = _mm_aesenc_si128    (m, key_schedule[ 1]);
    m = _mm_aesenc_si128    (m, key_schedule[ 2]);
    m = _mm_aesenc_si128    (m, key_schedule[ 3]);
    m = _mm_aesenc_si128    (m, key_schedule[ 4]);
    m = _mm_aesenc_si128    (m, key_schedule[ 5]);
    m = _mm_aesenc_si128    (m, key_schedule[ 6]);
    m = _mm_aesenc_si128    (m, key_schedule[ 7]);
    m = _mm_aesenc_si128    (m, key_schedule[ 8]);
    m = _mm_aesenc_si128    (m, key_schedule[ 9]);
    m = _mm_aesenclast_si128(m, key_schedule[10]);
    _mm_storeu_si128((__m128i *) cipherText, m);
  }
}


/*! \brief Two AES-ECB encryption with two distinct keys  for the GGM PRF */
void aes128_dbl_enc(u64 *plainText, u64 *cipherText, u64 nbblocks, __m128i key_schedule0[11], __m128i key_schedule1[11]){
  __m128i m, ms;

  for(u64 i=0;i<nbblocks;i++,plainText+=2,cipherText+=4) {
    m = _mm_loadu_si128((__m128i *) plainText);
    ms=m;
    m = _mm_xor_si128       (m, key_schedule0[ 0]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 1]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 2]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 3]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 4]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 5]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 6]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 7]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 8]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 9]);
    m = _mm_aesenclast_si128(m, key_schedule0[10]);

    m = _mm_xor_si128       (m, ms);
    _mm_storeu_si128((__m128i *) cipherText, m);
    
    m = ms; //_mm_loadu_si128((__m128i *) plainText);
    m = _mm_xor_si128       (m, key_schedule1[ 0]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 1]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 2]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 3]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 4]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 5]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 6]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 7]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 8]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 9]);
    m = _mm_aesenclast_si128(m, key_schedule1[10]);

    m = _mm_xor_si128       (m, ms);
    _mm_storeu_si128((__m128i *) (cipherText+2), m);
  }
}

/*! \brief Two AES-ECB encryption with two distinct keys tweaked into a delta-preserving doubling for the GGM PRF */
void aes128_dbl_enc_delta_preserve(u64 *plainText, u64 *cipherText, u64 nbblocks, __m128i key_schedule0[11], __m128i key_schedule1[11]){
  __m128i m, m1, m2, ms;

  for(u64 i=0;i<nbblocks;i++,plainText+=2,cipherText+=4) {
    m = _mm_loadu_si128((__m128i *) plainText);
    ms=m;
    m = _mm_xor_si128       (m, key_schedule0[ 0]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 1]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 2]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 3]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 4]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 5]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 6]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 7]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 8]);
    m = _mm_aesenc_si128    (m, key_schedule0[ 9]);
    m = _mm_aesenclast_si128(m, key_schedule0[10]);
    m1=m;
    
    m = ms; //_mm_loadu_si128((__m128i *) plainText);
    m = _mm_xor_si128       (m, key_schedule1[ 0]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 1]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 2]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 3]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 4]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 5]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 6]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 7]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 8]);
    m = _mm_aesenc_si128    (m, key_schedule1[ 9]);
    m2 = _mm_aesenclast_si128(m, key_schedule1[10]);

    
    m = _mm_xor_si128       (m1, m2);
    _mm_storeu_si128((__m128i *) cipherText, m);
    m = _mm_xor_si128       (m, ms);
    _mm_storeu_si128((__m128i *) (cipherText+2), m);
  }
}


/*! \brief Multiple AES-ECB encryptions with a single key and rotated inputs for final expansion */
void aes128_expand(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, __m128i key_schedule[11])
{
  __m128i m, m0, mask;

  // Product of two cyclic permutations to get high order 7*9=63 max replications
  // If need, resplit 9 into 4+5 and get up to 7*4*5=140 
  mask=_mm_setr_epi8(1,2,3,4,5,6,0,8,9,10,11,12,13,14,15,7);
  for(u64 loop=0;loop<nbblocks;loop++,plainText+=2){
    m = _mm_loadu_si128((__m128i *) plainText);
    m0=m;
    for(u64 i=0;i<blocksize;i++,cipherText+=2){
      m=m0;
      //      m0 = _mm_shuffle_epi8(m0, mask);
      m = _mm_xor_si128       (m, key_schedule[ 0]);
      m = _mm_aesenc_si128    (m, key_schedule[ 1]);
      m = _mm_aesenc_si128    (m, key_schedule[ 2]);
      m = _mm_aesenc_si128    (m, key_schedule[ 3]);
      m = _mm_aesenc_si128    (m, key_schedule[ 4]);
      m = _mm_aesenc_si128    (m, key_schedule[ 5]);
      m = _mm_aesenc_si128    (m, key_schedule[ 6]);
      m = _mm_aesenc_si128    (m, key_schedule[ 7]);
      m = _mm_aesenc_si128    (m, key_schedule[ 8]);
      m = _mm_aesenc_si128    (m, key_schedule[ 9]);
      m = _mm_aesenclast_si128(m, key_schedule[10]);
      m = _mm_xor_si128       (m, m0);
      m0 = _mm_shuffle_epi8(m0, mask);
      _mm_storeu_si128((__m128i *) (cipherText), m);
    }
  }
}

/*! \brief Multiple AES-ECB encryptions with a single key and rotated inputs for final expansion. First output block is a copy an the input. */
void aes128_expand_1cpy(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, __m128i key_schedule[11])
{
  __m128i m, m0, mask;

  mask=_mm_setr_epi8(1,2,3,4,5,6,0,8,9,10,11,12,13,14,15,7);
  for(u64 loop=0;loop<nbblocks;loop++,plainText+=2){
    m = _mm_loadu_si128((__m128i *) plainText);
    m0=m;
    _mm_storeu_si128((__m128i *) (cipherText), m);
    cipherText+=2;
   for(u64 i=1;i<blocksize;i++,cipherText+=2){
      m=m0;
      //      m0 = _mm_shuffle_epi8(m0, mask);
      m = _mm_xor_si128       (m, key_schedule[ 0]);
      m = _mm_aesenc_si128    (m, key_schedule[ 1]);
      m = _mm_aesenc_si128    (m, key_schedule[ 2]);
      m = _mm_aesenc_si128    (m, key_schedule[ 3]);
      m = _mm_aesenc_si128    (m, key_schedule[ 4]);
      m = _mm_aesenc_si128    (m, key_schedule[ 5]);
      m = _mm_aesenc_si128    (m, key_schedule[ 6]);
      m = _mm_aesenc_si128    (m, key_schedule[ 7]);
      m = _mm_aesenc_si128    (m, key_schedule[ 8]);
      m = _mm_aesenc_si128    (m, key_schedule[ 9]);
      m = _mm_aesenclast_si128(m, key_schedule[10]);
      m = _mm_xor_si128       (m, m0);
      m0 = _mm_shuffle_epi8(m0, mask);
      _mm_storeu_si128((__m128i *) (cipherText), m);

    }
  }
}
#else

void aes128_expand_key(aessubkeytype key_schedule[11], u64 enc_key[2]){
  aes_key_setup((const BYTE *) enc_key, (WORD *) key_schedule, 128);
}

void aes128_enc(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule[11]){
  for(u64 i=0;i<nbblocks;i++,plainText+=2,cipherText+=2) {
    aes_encrypt((const BYTE *) plainText, ( BYTE *) cipherText, (const WORD *) key_schedule, 128);
  }
}


void aes128_dbl_enc(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule0[11], aessubkeytype key_schedule1[11]){
  for(u64 i=0;i<nbblocks;i++,plainText+=2,cipherText+=4) {
    aes_encrypt((const BYTE *) plainText, ( BYTE *) cipherText, (const WORD *) key_schedule0, 128);
    aes_encrypt((const BYTE *) plainText, ( BYTE *) (cipherText+2), (const WORD *) key_schedule1, 128);
    cipherText[0]^=plainText[0];
    cipherText[1]^=plainText[1];
    cipherText[2]^=plainText[0];
    cipherText[3]^=plainText[1];
  }
}

void aes128_dbl_enc_delta_preserve(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule0[11], aessubkeytype key_schedule1[11]){
  for(u64 i=0;i<nbblocks;i++,plainText+=2,cipherText+=4) {
    aes_encrypt((const BYTE *) plainText, ( BYTE *) cipherText, (const WORD *) key_schedule0, 128);
    aes_encrypt((const BYTE *) plainText, ( BYTE *) (cipherText+2), (const WORD *) key_schedule1, 128);
    cipherText[0]^=cipherText[2];
    cipherText[1]^=cipherText[3];
    cipherText[2]=plainText[0]^cipherText[0];
    cipherText[3]=plainText[1]^cipherText[1];
  }
}

void aes128_expand(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, aessubkeytype key_schedule[11])
{
  u64 plainrotated[2], tmp;

  for(u64 loop=0;loop<nbblocks;loop++,plainText+=2){
    plainrotated[0]=plainText[0];
    plainrotated[1]=plainText[1];
    for(u64 i=0;i<blocksize;i++,cipherText+=2){
      aes_encrypt((const BYTE *) plainrotated, ( BYTE *) cipherText, (const WORD *) key_schedule, 128);
      tmp=plainrotated[0];
      plainrotated[0]>>=8;
      plainrotated[0]^=(plainrotated[1]<<56);
      plainrotated[1]>>=8;
      plainrotated[1]^=(tmp<<56);
    }
  }
}

void aes128_expand_1cpy(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, aessubkeytype key_schedule[11])
{
  u64 plainrotated[2], tmp;

  for(u64 loop=0;loop<nbblocks;loop++,plainText+=2){
    plainrotated[0]=plainText[0];
    plainrotated[1]=plainText[1];
    cipherText[0]=plainText[0];
    cipherText[1]=plainText[1];
    cipherText+=2;
    for(u64 i=1;i<blocksize;i++,cipherText+=2){
      aes_encrypt((const BYTE *) plainrotated, ( BYTE *) cipherText, (const WORD *) key_schedule, 128);
      tmp=plainrotated[0];
      plainrotated[0]>>=8;
      plainrotated[0]^=(plainrotated[1]<<56);
      plainrotated[1]>>=8;
      plainrotated[1]^=(tmp<<56);
    }
  }
}
#endif

/*! \brief SHA256 variant for the GGM PRF */
void sha256_dbl_enc(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule0[11], aessubkeytype key_schedule1[11]){
  SHA256_CTX hash_ctx0,  hash_ctx;
  
  SHA256Init(&hash_ctx0);
  SHA256Update(&hash_ctx0, (uchar *)key_schedule0, 16);

  for(u64 i=0;i<nbblocks;i++,plainText+=2,cipherText+=4) {
    hash_ctx=hash_ctx0;
    SHA256Update(&hash_ctx, (uchar *)(&i), sizeof(u64));
    SHA256Update(&hash_ctx, (uchar *)plainText, 16);
    SHA256Final((uchar *)cipherText, &hash_ctx);
  }
}

/*! \brief SHA256 variant (XOR preserving) for the GGM PRF */
void sha256_dbl_enc_delta_preserve(u64 *plainText, u64 *cipherText, u64 nbblocks, aessubkeytype key_schedule0[11], aessubkeytype key_schedule1[11]){
  SHA256_CTX hash_ctx0,  hash_ctx;
  
  SHA256Init(&hash_ctx0);
  SHA256Update(&hash_ctx0, (uchar *)key_schedule0, 16);

  for(u64 i=0;i<nbblocks;i++,plainText+=2,cipherText+=4) {
    hash_ctx=hash_ctx0;
    SHA256Update(&hash_ctx, (uchar *)(&i), sizeof(u64));
    SHA256Update(&hash_ctx, (uchar *)plainText, 16);
    SHA256Final((uchar *)cipherText, &hash_ctx);
    cipherText[0]^=cipherText[2];
    cipherText[1]^=cipherText[3];
    cipherText[2]=cipherText[0]^plainText[0];
    cipherText[3]=cipherText[1]^plainText[1];
  }
}

/*! \brief SHA256 variant for the GGM PRF expansion of leaves*/
void sha256_expand(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, aessubkeytype key_schedule[11])
{
  SHA256_CTX hash_ctx0,  hash_ctx1, hash_ctx;
  
  SHA256Init(&hash_ctx0);
  SHA256Update(&hash_ctx0, (uchar *)key_schedule, 16);

  for(u64 loop=0;loop<nbblocks;loop++,plainText+=2){
    hash_ctx1=hash_ctx0;
    SHA256Update(&hash_ctx1, (uchar *)plainText, 16);
    if (blocksize&1) {
      hash_ctx=hash_ctx1;
      SHA256Final((uchar *)cipherText, &hash_ctx);
      cipherText+=2;
    }
    for(u64 i=blocksize&1;i<blocksize;i+=2,cipherText+=4){
      hash_ctx=hash_ctx1;
      SHA256Update(&hash_ctx, (uchar *)(&i), sizeof(u64));
      SHA256Final((uchar *)cipherText, &hash_ctx);
    }
  }
}

/*! \brief SHA256 variant for the GGM PRF expansion of leaves (XOR preserving)*/
void sha256_expand_1cpy(u64 *plainText, u64 *cipherText, u64 blocksize, u64 nbblocks, aessubkeytype key_schedule[11])
{
  SHA256_CTX hash_ctx0,  hash_ctx1, hash_ctx;
  
  SHA256Init(&hash_ctx0);
  SHA256Update(&hash_ctx0, (uchar *)key_schedule, 16);
  blocksize--;
  
  for(u64 loop=0;loop<nbblocks;loop++,plainText+=2){
    cipherText[0]=plainText[0];
    cipherText[1]=plainText[1];
    cipherText+=2;
    hash_ctx1=hash_ctx0;
    SHA256Update(&hash_ctx1, (uchar *)plainText, 16);
    if (blocksize&1) {
      hash_ctx=hash_ctx1;
      SHA256Final((uchar *)cipherText, &hash_ctx);
      cipherText+=2;
    }
    for(u64 i=blocksize&1;i<blocksize;i+=2,cipherText+=4){
      hash_ctx=hash_ctx1;
      SHA256Update(&hash_ctx, (uchar *)(&i), sizeof(u64));
      SHA256Final((uchar *)cipherText, &hash_ctx);
    }
  }
}



/*! \brief Timing init 

All timing functions are taken from the SDiTH public implementation*/
void btimer_init(btimer_t* timer) {
    if(timer != NULL) {
        timer->counter = 0;
        timer->nb_milliseconds = 0.;
        timer->nb_cycles = 0;
    }
}
void btimer_count(btimer_t *timer) {
    if(timer != NULL)
        timer->counter++;
}
void btimer_start(btimer_t *timer) {
    if(timer != NULL) {
        gettimeofday(&timer->start, NULL);
      #ifdef BENCHMARK_CYCLES
        timer->cstart = __rdtscp(&timer->garbage);
      #endif
    }
}

double btimer_diff(btimer_t *timer) {
    return ( (timer->stop.tv_sec - timer->start.tv_sec) * 1000000 + (timer->stop.tv_usec - timer->start.tv_usec) )/1000.;
}
uint64_t btimer_diff_cycles(btimer_t *timer) {
    return (timer->cstop - timer->cstart);
}
void btimer_end(btimer_t *timer) {
    if(timer != NULL) {
        gettimeofday(&timer->stop, NULL);
      #ifdef BENCHMARK_CYCLES
        timer->cstop = __rdtscp(&timer->garbage);
      #endif
        timer->nb_milliseconds += btimer_diff(timer);
      #ifdef BENCHMARK_CYCLES
        timer->nb_cycles += btimer_diff_cycles(timer);
      #endif
    }
}
double btimer_get(btimer_t *timer) {
    return timer->nb_milliseconds/timer->counter;
}
double btimer_get_cycles(btimer_t *timer) {
    return (double)timer->nb_cycles/timer->counter;
}





