#include <stdlib.h>
#include <string.h>

#include "assertions.h"
#include "param.h"
#include "rng.h"
#include "types.h"

typedef struct treeprg_sha3_context_struct {
  salt_t salt;
} treeprg_sha3_context_t;

typedef struct seed_expand_vec {
  salt_t salt;
  uint16_t iteration;
  uint16_t node_idx;
  seed_t parent_seed;
} seed_expand_vec_t;

#ifdef TEST_AES
__m128i aes_128_key_expansion(__m128i key, __m128i keygened){
	keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	return _mm_xor_si128(key, keygened);
}

//public API
EXPORT void aes128_load_key(__m128i key_schedule[11], uint64_t enc_key[2]){
    key_schedule[0] = _mm_set_epi64x(enc_key[0], enc_key[1]);
	key_schedule[1]  = AES_128_key_exp(key_schedule[0], 0x01);
	key_schedule[2]  = AES_128_key_exp(key_schedule[1], 0x02);
	key_schedule[3]  = AES_128_key_exp(key_schedule[2], 0x04);
	key_schedule[4]  = AES_128_key_exp(key_schedule[3], 0x08);
	key_schedule[5]  = AES_128_key_exp(key_schedule[4], 0x10);
	key_schedule[6]  = AES_128_key_exp(key_schedule[5], 0x20);
	key_schedule[7]  = AES_128_key_exp(key_schedule[6], 0x40);
	key_schedule[8]  = AES_128_key_exp(key_schedule[7], 0x80);
	key_schedule[9]  = AES_128_key_exp(key_schedule[8], 0x1B);
	key_schedule[10] = AES_128_key_exp(key_schedule[9], 0x36);
}

void aes128_enc(uint64_t plainText[2], uint64_t cipherText[2], __m128i key_schedule[11]){
    __m128i m = _mm_loadu_si128((__m128i *) plainText);

    DO_ENC_BLOCK(m, key_schedule);

    _mm_storeu_si128((__m128i *) cipherText, m);
}


/** @brief takes as input n/2 seeds and produces n seeds using AES:
 * out[i]=F_{first_tweak+i}(in[i/2]) */
EXPORT void sdith_tree_aes_seed_expand(TREE_PRG_CTX *key, void *out,
                                       const void *in,
                                       const uint16_t first_tweak,
                                       const uint16_t iteration, uint64_t n,
                                       __m128i key_schedule_k0[11], 
                                       __m128i key_schedule_k1[11]) {
  ASSERT_DRAMATICALLY(n % 2 == 0, "n must be even for this function");
  treeprg_sha3_context_t *context = (treeprg_sha3_context_t *)key;
  const seed_t *in_seeds = (seed_t *)in;
  seed_t(*out_seeds)[2] = (seed_t(*)[2])out;

  uint64_t in_aes[2] = {0,0};
  uint64_t out_aes[2] = {0,0};

  for (uint64_t i = 0; i < n / 2; ++i) {
    // out[i, i+1] = aes(in[i], K_0), aes(in[i], K_1)
    memcpy(in_aes, in_seeds + i, PARAM_seed_size);
    aes128_enc(in_aes, out_aes, key_schedule_k0);
    out_aes[0] ^= in_aes[0];
    out_aes[1] ^= in_aes[1];
    memcpy(out_seeds[i][0], out_aes, 16);
    aes128_enc(in_aes, out_aes, key_schedule_k1);
    out_aes[0] ^= in_aes[0];
    out_aes[1] ^= in_aes[1];
    memcpy(out_seeds[i][1], out_aes, 16);
  }
}
#endif

/** @brief takes as input n/2 seeds and produces n seeds:
 * out[i]=F_{first_tweak+i}(in[i/2]) */
EXPORT void sdith_tree_prg_seed_expand(TREE_PRG_CTX *key, void *out,
                                       const void *in,
                                       const uint16_t first_tweak,
                                       const uint16_t iteration, uint64_t n) {
  ASSERT_DRAMATICALLY(n % 2 == 0, "n must be even for this function");
  treeprg_sha3_context_t *context = (treeprg_sha3_context_t *)key;
  const seed_t *in_seeds = (seed_t *)in;
  seed_t(*out_seeds)[2] = (seed_t(*)[2])out;

  // the size of the current layer is also its starting index.
  uint32_t tweak = first_tweak;

  if (n < 8) {
    seed_expand_vec_t vec;
    memcpy(vec.salt, context->salt, PARAM_salt_size);
    vec.iteration = iteration;

    for (uint64_t i = 0; i < n / 2; ++i) {
      // out[i, i+1] = sha3(salt||iteration||node_idx||in[i])
      vec.node_idx = tweak;
      memcpy(vec.parent_seed, in_seeds + i, PARAM_seed_size);
      sdith_hash(HASH_TREE, out_seeds[i], &vec, sizeof(vec));
      ++tweak;
    }
    return;
  }

  seed_expand_vec_t vec4[4];
  void *src4[4] = {&vec4[0], &vec4[1], &vec4[2], &vec4[3]};
  for (uint64_t i = 0; i < 4; ++i) {
    memcpy(vec4[i].salt, context->salt, PARAM_salt_size);
    vec4[i].iteration = iteration;
  }
  for (uint64_t i = 0; i < n / 2; i += 4) {
    vec4[0].node_idx = tweak;
    memcpy(vec4[0].parent_seed, in_seeds + i, PARAM_seed_size);
    vec4[1].node_idx = tweak + 1;
    memcpy(vec4[1].parent_seed, in_seeds + i + 1, PARAM_seed_size);
    vec4[2].node_idx = tweak + 2;
    memcpy(vec4[2].parent_seed, in_seeds + i + 2, PARAM_seed_size);
    vec4[3].node_idx = tweak + 3;
    memcpy(vec4[3].parent_seed, in_seeds + i + 3, PARAM_seed_size);

    void *dst4[4] = {out_seeds[i], out_seeds[i + 1], out_seeds[i + 2],
                     out_seeds[i + 3]};
    sdith_hash4(HASH_TREE, dst4, src4, sizeof(seed_expand_vec_t));
  }
}

/** @brief takes a tree_prg context */
EXPORT TREE_PRG_CTX *sdith_create_tree_prg_ctx(void const *const root_salt) {
  treeprg_sha3_context_t *context =
      (treeprg_sha3_context_t *)malloc(sizeof(treeprg_sha3_context_t));
  memcpy(context->salt, root_salt, PARAM_salt_size);
  return (TREE_PRG_CTX *)context;
}

/** @brief free a tree_prg context */
EXPORT void sdith_free_tree_prg_ctx(TREE_PRG_CTX *key) { free(key); }

/** @brief expand the leaf into seed and rho */
EXPORT void sdith_tree_prg_leaf_expand(TREE_PRG_CTX *key, seed_t in_seed,
                                       seed_t out_seed,
                                       uint8_t out_rho[PARAM_rho_size]) {
  treeprg_sha3_context_t *context = (treeprg_sha3_context_t *)key;
  uint8_t in[PARAM_seed_size + PARAM_salt_size];
  memcpy(in, in_seed, PARAM_seed_size);
  memcpy(&in[PARAM_seed_size], context->salt, PARAM_salt_size);

  uint8_t out[PARAM_seed_size + PARAM_rho_size];
  sdith_hash(HASH_TREE, out, in, sizeof(in));
  memcpy(out_seed, out, PARAM_seed_size);
  memcpy(out_rho, &out[PARAM_seed_size], PARAM_rho_size);
}

/** @brief expand the leaf into seed and rho x4 */
EXPORT void sdith_tree_prg_leaf_expand4(TREE_PRG_CTX *key, seed_t *in_seed,
                                        seed_t *out_seed, uint8_t **out_rho) {
  treeprg_sha3_context_t *context = (treeprg_sha3_context_t *)key;
  uint8_t in4[4][PARAM_seed_size + PARAM_salt_size];
  for (uint64_t i = 0; i < 4; ++i) {
    memcpy(in4[i], in_seed[i], PARAM_seed_size);
    memcpy(&in4[i][PARAM_seed_size], context->salt, PARAM_salt_size);
  }
  void *inVec[4] = {&in4[0], &in4[1], &in4[2], &in4[3]};
  uint8_t out4[4][PARAM_seed_size + PARAM_rho_size];
  void *outVec[4] = {&out4[0], &out4[1], &out4[2], &out4[3]};

  sdith_hash4(HASH_TREE, outVec, inVec, sizeof(in4[0]));

  for (uint64_t i = 0; i < 4; ++i) {
    memcpy(out_seed[i], out4[i], PARAM_seed_size);
    memcpy(out_rho[i], &out4[i][PARAM_seed_size], PARAM_rho_size);
  }
}
