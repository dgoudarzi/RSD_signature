#include "pprf_sig.h"

#define B_TREE_GEN 3

static __m128i aes_128_key_expansion(__m128i key, __m128i keygened){
	keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	return _mm_xor_si128(key, keygened);
}

//public API
void aes128_load_key(__m128i key_schedule[11], u64 enc_key[2]){
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

void aes128_enc(u64 plainText[2], u64 cipherText[2], __m128i key_schedule[11]){
    __m128i m = _mm_loadu_si128((__m128i *) plainText);

    DO_ENC_BLOCK(m, key_schedule);

    _mm_storeu_si128((__m128i *) cipherText, m);
}

// Randomness functions

void sample_random (u64* dest, size_t nb_bits, __m128i key_schedule[11]) {
    size_t nb_iter = (nb_bits/128);
    size_t reminder = nb_bits - (nb_iter*128);
    unsigned long long r0, r1;
    _rdrand64_step(&r0);
    r1 = 0;
    unsigned long long m0, m1;
    _rdrand64_step(&m0);
    _rdrand64_step(&m1);
    u64 rd_in[2] = {r0, r1};
    u64 tmp[2] = {m0, m1};

    aes128_enc(rd_in,  &dest[0], key_schedule);
    dest[0] ^= tmp[0];
    dest[1] ^= tmp[1];
    rd_in[1]++;
    for (size_t i=1; i<nb_iter; i++) {
        aes128_enc(rd_in,  &dest[i*2], key_schedule);
        dest[i*2] ^= dest[(i-1)*2];
        dest[i*2+1] ^= dest[(i-1)*2+1];
        rd_in[1]++;
    }
    if (reminder != 0) {
        aes128_enc(rd_in, tmp, key_schedule);
        if (reminder < 64) {
            dest[nb_iter*2] = (u64)(tmp[0] ^ dest[(nb_iter-1)*2]) >> (64 - reminder);
        }
        else if (reminder < 128) {
            dest[nb_iter*2] = tmp[0] ^ dest[(nb_iter-1)*2];
            dest[(nb_iter*2)+1] = (u64)(tmp[0] ^ dest[((nb_iter-1)*2)+1]) >> (128 - reminder);
        }
    }
}

void sample_random_from_seed (u64 seed[2], u64* dest, size_t nb_bits, __m128i key_schedule[11]) {
    size_t nb_iter = (nb_bits/128);
    size_t reminder = nb_bits - (nb_iter*128);
    unsigned long long r0, r1;
    r0 = seed[0];
    r1 = 0;
    u64 rd_in[2] = {r0, r1};
    
    aes128_enc(rd_in,  &dest[0], key_schedule);
    dest[0] ^= seed[0];
    dest[1] ^= seed[1];
    rd_in[1]++;
    for (size_t i=1; i<nb_iter; i++) {
        aes128_enc(rd_in,  &dest[i*2], key_schedule);
        dest[i*2] ^= dest[(i-1)*2];
        dest[i*2+1] ^= dest[(i-1)*2+1];
        rd_in[1]++;
    }
    u64 tmp[2];
    if (reminder != 0) {
        aes128_enc(rd_in, tmp, key_schedule);
        if (reminder < 64) {
            dest[nb_iter*2] = (u64)(tmp[0] ^ dest[(nb_iter-1)*2]) >> (64 - reminder);
        }
        else if (reminder < 128) {
            dest[nb_iter*2] = tmp[0] ^ dest[(nb_iter-1)*2];
            dest[(nb_iter*2)+1] = (u64)(tmp[0] ^ dest[((nb_iter-1)*2)+1]) >> (128 - reminder);
        }
    }
}

// Utils
#define extract_and_power_two_of(x, y) (1<<(((x)>>(4*(y))) & 0xF))
#define extract_four_bits(x, y) (((x) >>(4*(y))) & 0xF)
#define extract_sixteen_bits(x,y) (((x) >> (16*(y))) & 0xFFFF)

uint16_t rotl16 (uint16_t value, unsigned int count) {
    const unsigned int mask = CHAR_BIT * sizeof(value) - 1;
    count &= mask;
    return (value << count) | (value >> (-count & mask));
}

void expand(u64 src[word_x], u64 dst[word_K]) {
    for (size_t i=0; i < word_x-1; i++) {
        for (size_t pos=0; pos<16; pos++) {
            size_t dst_pos = pos/4;
            size_t shift_pos = pos % 4;
            dst[i*4 + dst_pos] ^= (u64) extract_and_power_two_of(src[i], pos) << (16*shift_pos);
        }
    }
    for (size_t pos=0; pos<13; pos++) {
        size_t dst_pos = pos/4;
        size_t shift_pos = pos %4;
        dst[(word_x-1)*4 + dst_pos] ^= (u64) extract_and_power_two_of(src[(word_x-1)], pos) << (16*shift_pos);
    }
}

// a optimiser avec des instructions vectorielles (last truc a faire)
void mul_mat_vec(u64 dest[word_k], u64 mat[word_H], u64 vec[word_K]) {
    for (size_t i=0; i<word_k; i++) {
        dest[i] = 0;
        for (size_t j=0; j<word_K; j++) {
            dest[i] ^= vec[j] & mat[(i*word_K) + j];
        }
    }
}

size_t log2_d(size_t in) {
    size_t r = 0; // r will be lg(v)
    while (in >>= 1){
        r++;    
    }
    return r;
}

#ifdef PPRF_SHA

#include "sha3/KeccakHash.h"
#include "sha3/KeccakHashtimes4.h"

void SHA_sdith(u64 dest[2], u64 parent[2], size_t index_leaf, size_t index_e, u64 salt[word_lambda]) {
    unsigned char buff_h1[SHA256_DIGEST_LENGTH];
    Keccak_HashInstance *h1 =
      (Keccak_HashInstance *)malloc(sizeof(Keccak_HashInstance));
    Keccak_HashInitialize_SHA3_256(h1);
    Keccak_HashUpdate(h1, (uint8_t  *) parent, 2*64/8);
    Keccak_HashUpdate(h1, (uint8_t  *) &index_leaf, 1);
    Keccak_HashUpdate(h1, (uint8_t  *) &index_e, 1);
    Keccak_HashUpdate(h1, (uint8_t  *) salt, lambda/8);
    Keccak_HashFinal(h1, buff_h1);
    memcpy(dest, buff_h1, 2*64/8);
}

 void pprf_sha(u64 dest[2], u64 path[D][2], size_t index_leaf, size_t e, u64 K_0[word_lambda], u64 K_1[word_lambda]) {
    if (index_leaf == 0) {
        for (size_t i=1; i<D; i++) {
            SHA_sdith(path[i], path[i-1], index_leaf, e, K_0);
        }
        SHA_sdith(dest, path[D-1], index_leaf, e, K_0);
    }
    else if (index_leaf & 1) {
        SHA_sdith(dest, path[D-1], index_leaf, e, K_1);
    }
    else {
        size_t lg_diff = log2_d(index_leaf ^ (index_leaf-1));
        size_t restart_floor = D-lg_diff;
        for (size_t i=restart_floor; i<D; i++) {
            if ((index_leaf>>(D-i)) & 1) {
                SHA_sdith(path[i], path[i-1], index_leaf, e, K_1);
            }
            else {
                SHA_sdith(path[i], path[i-1], index_leaf, e, K_0);
            }
        }
        SHA_sdith(dest, path[D-1], index_leaf, e, K_0);
    }
 }

 void compute_copath_sha(u64 dest[D][2], u64 in[2], size_t index_leaf, size_t e, u64 K_0[word_lambda], u64 K_1[word_lambda]) {
    u64 tmp[2];
    for (size_t d=1; d<D; d++) {
        if ((index_leaf>>(D-d)) & 1) {
            SHA_sdith(tmp, in, index_leaf, e, K_1);
            SHA_sdith(dest[d-1], in, index_leaf, e, K_0);
        }
        else{
            SHA_sdith(tmp, in, index_leaf, e, K_0);
            SHA_sdith(dest[d-1], in, index_leaf, e, K_1);
        }
        in[0] = tmp[0];
        in[1] = tmp[1];    
    }
    if (index_leaf & 1) {
        SHA_sdith(dest[D-1], in, index_leaf, e, K_0);
    }   
    else {
        SHA_sdith(dest[D-1], in, index_leaf, e, K_1);
    }
 }

void recover_seed_from_copath_sha(u64 dest[2], u64 copath[D][2], size_t index_ie, size_t index_leaf, size_t e, u64 K_0[word_lambda], u64 K_1[word_lambda]) {
    size_t start_floor = 0;
    u64 path[D][2]; 
    for (size_t d = 0; d<D; d++) {
        if ((index_leaf>>d) == ((index_ie>>d) ^ 1))  {
            start_floor = D-d; 
            path[start_floor-1][0] = copath[start_floor-1][0];
            path[start_floor-1][1] = copath[start_floor-1][1];
            for (size_t i=start_floor; i<D; i++) {
                if ((index_leaf>>(D-i-1)) & 1) {
                    SHA_sdith(path[i], path[i-1], index_leaf, e, K_1);
                }
                else {
                    SHA_sdith(path[i], path[i-1], index_leaf, e, K_0);
                }
            }
        }
    }
    dest[0] = path[D-1][0];
    dest[1] = path[D-1][1];
}

#else

void pprf(u64 dest[2], u64 path[D][2], size_t index_leaf, __m128i K_0[11], __m128i K_1[11]) {
    if (index_leaf == 0) {
        for (size_t i=1; i<D; i++) {
            aes128_enc(path[i-1], path[i], K_0);
            path[i][0] ^= path[i-1][0];
            path[i][1] ^= path[i-1][1];
        }
        aes128_enc(path[D-1], dest, K_0);
        dest[0] ^= path[D-1][0];
        dest[1] ^= path[D-1][1];
    }
    else if (index_leaf & 1) {
        aes128_enc(path[D-1], dest, K_1);
        dest[0] ^= path[D-1][0];
        dest[1] ^= path[D-1][1];
    }
    else {
        size_t lg_diff = log2_d(index_leaf ^ (index_leaf-1));
        size_t restart_floor = D-lg_diff;
        for (size_t i=restart_floor; i<D; i++) {
            if ((index_leaf>>(D-i)) & 1) {
                aes128_enc(path[i-1], path[i], K_1);
            }
            else {
                aes128_enc(path[i-1], path[i], K_0);
            }
            path[i][0] ^= path[i-1][0];
            path[i][1] ^= path[i-1][1];
        }
        aes128_enc(path[D-1], dest, K_0);
        dest[0] ^= path[D-1][0];
        dest[1] ^= path[D-1][1];
    }
 }

void compute_copath(u64 dest[D][2], u64 in[2], size_t index_leaf, __m128i K_0[11], __m128i K_1[11]) {
    u64 tmp[2];
    for (size_t d=1; d<D; d++) {
        if ((index_leaf>>(D-d)) & 1) {
            aes128_enc(in, tmp, K_1);
            aes128_enc(in, dest[d-1], K_0);
        }
        else{
            aes128_enc(in, tmp, K_0);
            aes128_enc(in, dest[d-1], K_1);
        }
        tmp[0] ^= in[0];
        tmp[1] ^= in[1];
        dest[d-1][0] ^= in[0];
        dest[d-1][1] ^= in[1];
        in[0] = tmp[0];
        in[1] = tmp[1];    
    }
    if (index_leaf & 1) {
        aes128_enc(in, dest[D-1], K_0);
    }   
    else {
        aes128_enc(in, dest[D-1], K_1);
    }
    dest[D-1][0] ^= in[0];
    dest[D-1][1] ^= in[1];

 }

void recover_seed_from_copath(u64 dest[2], u64 copath[D][2], size_t index_ie, size_t index_leaf, __m128i K_0[11], __m128i K_1[11]) {
    size_t start_floor = 0;
    u64 path[D][2]; 
    for (size_t d = 0; d<D; d++) {
        if ((index_leaf>>d) == ((index_ie>>d) ^ 1))  {
            start_floor = D-d; 
            path[start_floor-1][0] = copath[start_floor-1][0];
            path[start_floor-1][1] = copath[start_floor-1][1];
            for (size_t i=start_floor; i<D; i++) {
                if ((index_leaf>>(D-i-1)) & 1) {
                    aes128_enc(path[i-1], path[i], K_1);
                }
                else {
                    aes128_enc(path[i-1], path[i], K_0);
                }
                path[i][0] ^= path[i-1][0];
                path[i][1] ^= path[i-1][1];
            }
        }
    }
    dest[0] = path[D-1][0];
    dest[1] = path[D-1][1];
}
#endif

int crypto_sign_keypair (u64 sk_x[word_x], 
            u64 seed[word_lambda], 
            u64 K_0[word_lambda], 
            u64 K_1[word_lambda], 
            u64 H[word_H], 
            u64 y[word_k], 
            __m128i key_schedule_k0[11], 
            __m128i key_schedule_k1[11]
    ) {
    // Sample (K_0, K_1) <-_r {0,1}^{lambda}x{0,1}^{lambda}
    unsigned long long kl, kh;
	_rdrand64_step(&kh);
	_rdrand64_step(&kl);
    K_0[0] = kh;
    K_0[1] = kl;
    aes128_load_key(key_schedule_k0, K_0);
   	_rdrand64_step(&kh);
	_rdrand64_step(&kl);
    K_1[0] = kh;
    K_1[1] = kl;
    aes128_load_key(key_schedule_k1, K_1);
    // Set H <- PRG(seed) and y <- H * Expand(x)
    sample_random_from_seed(seed, H, (size_t) (big_K*small_k), key_schedule_k0);
    sample_random(sk_x, (size_t)(log_bs*w), key_schedule_k0);
    u64 expand_x[word_K];
    for (size_t k=0; k<word_K; k++) {
       expand_x[k] = 0;
    }
    expand(sk_x, expand_x);
    mul_mat_vec(y, H, expand_x);

    return 0;
}

// Signature algorith of the signature scheme 
int crypto_sign (u64 *m, 
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
    ) {
    //// Initialization
    unsigned int resultLen = 0;
    unsigned char buff_h1[SHA256_DIGEST_LENGTH];
    unsigned char buff_h2[SHA256_DIGEST_LENGTH];
    EVP_MD_CTX* h1 = EVP_MD_CTX_new();
    EVP_DigestInit_ex(h1, EVP_sha256(), NULL);
    EVP_MD_CTX* h2 = EVP_MD_CTX_new();
    EVP_DigestInit_ex(h2, EVP_sha256(), NULL);

    EVP_DigestUpdate(h1, m, m_len);
    EVP_DigestUpdate(h1, K_0, lambda/8);
    EVP_DigestUpdate(h1, K_1, lambda/8);
    
    EVP_DigestUpdate(h2, m, m_len);
    EVP_DigestUpdate(h2, K_0, lambda/8);
    EVP_DigestUpdate(h2, K_1, lambda/8);

    //// Phase 1
    u64 seed_e[tau][word_lambda];
    u64 X[tau][D][word_x], R[tau][D][word_x], U[tau][D][word_K];
    u64 x[tau][word_x], u[tau][word_K], r[tau][word_x];
    u64 r_n[tau][word_x];
    u64 ue[word_K];
    u64 tree[D][2];

    for (size_t e=0; e<tau; e++) {
        for (size_t d=0; d<D; d++) {
            for (size_t k=0; k<word_x; k++) {
                X[e][d][k] = 0;
                R[e][d][k] = 0;
            }
            for (size_t k=0; k<word_K; k++) {
                U[e][d][k] = 0;      
            }
        }
        for (size_t k=0; k<word_x; k++) {
            x[e][k] = 0;
            r[e][k] = 0;
            r_n[e][k] = 0;
        }
        for (size_t k=0; k<word_K; k++) {
            u[e][k] = 0;
        }
    }
    
    u64 seed_i[word_lambda];
    u64 sampled_data[((lambda + big_K + 2 * log_bs * w)/64)+ 8];
    u64 extracted_pos = 0;

    u64 log_com[tau][n_par-1][word_lambda];

    for (size_t e=0; e<tau; e++) {
        //// Initialization
        btimer_start(&timers_algos[B_TREE_GEN]);
        sample_random(seed_e[e], lambda, key_schedule_k0);
        for (size_t i=0; i<word_x; i++) {
            x[e][i] = sk_x[i];       
        }
        //// loop over parties
        tree[0][0] = seed_e[e][0];
        tree[0][1] = seed_e[e][1];
        btimer_end(&timers_algos[B_TREE_GEN]);
        for (size_t i=0; i<n_par-1; i++) {
            // Compute seed_i^e <- PPRF_salt(seed^e, i)
            btimer_start(&timers_algos[B_TREE_GEN]);
            #ifdef PPRF_SHA
            pprf_sha(seed_i, tree, i, e, K_0, K_1);
            #else
            pprf(seed_i, tree, i, key_schedule_k0, key_schedule_k1);
            #endif
            btimer_end(&timers_algos[B_TREE_GEN]);
            // (x_i^e, r_i^e, u_i^e, com_i^e) <- PRG(seed_i^e)
            sample_random_from_seed(seed_i, sampled_data, lambda + big_K + 2*(log_bs * w)+ 4*64, key_schedule_k0);
            // x_n^e <- x_n^e - x_i^e mod bs, u_n^e <- u_n^e xor u_i^e, and r^e <- r^e+r_i^e mod bs
            extracted_pos = 0;
            for (size_t k=0; k<word_x-1; k++) {
                for (size_t pos=0; pos<BS; pos++) {
                    x[e][k] = (x[e][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(x[e][k], pos) 
                                - extract_four_bits(sampled_data[k], pos)) % BS) << (4*pos));
                }
            }
            for (size_t pos=0; pos<13; pos++) {
                x[e][word_x-1] = (x[e][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(x[e][word_x-1], pos) 
                            - extract_four_bits(sampled_data[word_x-1], pos)) % BS) << (4*pos));
            }
            for (size_t k=0; k<word_x-1; k++) {
                for (size_t pos=0; pos<BS; pos++) {
                    r[e][k] = (r[e][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(r[e][k], pos) 
                                + extract_four_bits(sampled_data[word_x+k], pos)) % BS) << (4*pos));
                }
            }
            for (size_t pos=0; pos<13; pos++) {
                r[e][word_x-1] = (r[e][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(r[e][word_x-1], pos) 
                            + extract_four_bits(sampled_data[word_x+word_x-1], pos)) % BS) << (4*pos));
            }
            extracted_pos += 2*word_x;
            for (size_t k=0; k<word_K; k++) {
                u[e][k] ^= sampled_data[extracted_pos + k];
            }
            u[e][word_K-1] = u[e][word_K-1] & 0xFFFF;
            extracted_pos += word_K;
            // h(com_i^e)
            EVP_DigestUpdate(h1, &sampled_data[extracted_pos], lambda/8);
            memcpy(&log_com[e][i][0], &sampled_data[extracted_pos], word_lambda*8);
            
            // For all d <= D such that i[d] = 0, set
            for (size_t d=0; d<D; d++) {
                if (((i>>d) & 1) == 0) {
                    // X_d,0^e <- X_d,0^e + x_i^e mod bs
                    // R_d,0^e <- R_d,0^e + r_i^e mod bs
                    for (size_t k=0; k<word_x-1; k++) {
                        for (size_t pos=0; pos<BS; pos++) {
                            X[e][d][k] = (X[e][d][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[k], pos) 
                                        + extract_four_bits(X[e][d][k], pos)) % BS) << (4*pos));                           
                        }
                    }
                    for (size_t pos=0; pos<13; pos++) {
                            X[e][d][word_x-1] = (X[e][d][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x-1], pos) 
                                        + extract_four_bits(X[e][d][word_x-1], pos)) % BS) << (4*pos));
                    }
                    for (size_t k=0; k<word_x-1; k++) {
                        for (size_t pos=0; pos<BS; pos++) {
                            R[e][d][k] = (R[e][d][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x + k], pos) 
                                        + extract_four_bits(R[e][d][k], pos)) % BS) << (4*pos));
                           
                        }
                    }
                    for (size_t pos=0; pos<13; pos++) {
                            R[e][d][word_x-1] = (R[e][d][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x+word_x-1], pos) 
                                        + extract_four_bits(R[e][d][word_x-1], pos)) % BS) << (4*pos));
                    }
                    // U_d,0^e <- U_d,0^e xor u_i^e
                    for (size_t k=0; k<word_K; k++) {
                        U[e][d][k] ^= sampled_data[2*word_x + k];
                    }
                    U[e][d][word_K-1] = U[e][d][word_K-1] & 0xFFFF;
                }
            }
        }
        
        //// On node n:
        // Compute seed_n^e <- PPRF_salt(seed^e, n)
        btimer_start(&timers_algos[B_TREE_GEN]);
        #ifdef PPRF_SHA
        SHA_sdith(seed_i, tree[D-1], n_par-1, e, K_1);
        #else
        aes128_enc(tree[D-1], seed_i, key_schedule_k1);
        seed_i[0] ^= tree[D-1][0];
        seed_i[1] ^= tree[D-1][1];
        #endif
        btimer_end(&timers_algos[B_TREE_GEN]);
        // Compute r_n^e <- PRG(seed_n^e)
        sample_random_from_seed(seed_i, r_n[e], log_bs*w, key_schedule_k0);
        // r^e <- r^e + r_n^e mod bs, u^e <- Expand(r^e), and u_n^e<-u_n^e xor u^e
        for (size_t k=0; k<word_x-1; k++) {
            for (size_t pos=0; pos<BS; pos++) {
                r[e][k] = (r[e][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(r[e][k], pos) 
                            + extract_four_bits(r_n[e][k], pos)) % BS) << (4*pos));
            }
        }
        for (size_t pos=0; pos<13; pos++) {
            r[e][word_x-1] = (r[e][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(r[e][word_x-1], pos) 
                        + extract_four_bits(r_n[e][word_x-1], pos)) % BS) << (4*pos));
        }
        for (size_t k=0; k<word_K; k++) {
            ue[k] = 0;
        }
        expand(r[e], ue);
        for (size_t k=0; k<word_K; k++) {
            u[e][k] ^= ue[k];
        }
        // com_n^e <-H(state_n^e)
        EVP_MD_CTX* h = EVP_MD_CTX_new();
        EVP_DigestInit_ex(h, EVP_sha256(), NULL);
        EVP_DigestUpdate(h, &x[e][0], word_x*8);
        EVP_DigestUpdate(h, &u[e][0], word_K*8);
        EVP_DigestUpdate(h, seed_i, lambda/8);

        unsigned char buff_h[SHA256_DIGEST_LENGTH];
        EVP_DigestFinal_ex(h, buff_h, &resultLen);
        // Update h1 with com_n^e
        EVP_DigestUpdate(h1, buff_h, SHA256_DIGEST_LENGTH);
        memcpy(&sigma->com_n[e], buff_h, SHA256_DIGEST_LENGTH);
    }


    EVP_DigestFinal_ex(h1, buff_h1, &resultLen);
    memcpy(sigma->h1, buff_h1, SHA256_DIGEST_LENGTH);
    EVP_DigestUpdate(h2, buff_h1, SHA256_DIGEST_LENGTH);
    
    //// Phase 2
    // pi^e <- PRG_1(h1)
    uint8_t pi[w];
    uint8_t rand_pi[tau][w];
    sample_random_from_seed((uint64_t*)&buff_h1[0], (uint64_t*) &rand_pi[0][0], tau*w*8, key_schedule_k0);

    //// Phase 3
    // For each iteration e in tau
    u64 z[tau][word_x];
    u64 z_d[2][word_x];
    u64 y_d[2][word_k];
    u64 shift_u[word_K];
    size_t pos_pi;


    for (size_t e=0; e<tau; e++) {
        
        // generate random permutation using Fisher-Yates
        for (size_t i=0; i<w; i++) {
            pi[i] = i;
        }
        for (size_t i=w-1; i<w; --i) {
            size_t j = rand_pi[e][i] % (i+1);
            size_t tmp = pi[i];
            pi[i] = pi[j];
            pi[j] = tmp;
        }
        
        // z^e <- x - pi^e(r^e) mod bs  
        pos_pi = 0;
        for (size_t k=0; k<word_x-1; k++) {
            for (size_t pos=0; pos<BS; pos++) {
                z[e][k] = (z[e][k] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(sk_x[k], pos) 
                        - extract_four_bits(r[e][pi[pos_pi]/16], pi[pos_pi]%16)) % BS) << (4*pos));
                pos_pi++;
            }
        }
        for (size_t pos=0; pos<13; pos++) {
            z[e][word_x-1] = (z[e][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(sk_x[word_x-1], pos) 
                          - extract_four_bits(r[e][pi[pos_pi]/16], pi[pos_pi]%16)) % BS) << (4*pos));
            pos_pi++;
        }
        
        for (size_t d=0; d<D; d++) {
            // y_d0^e <- H * Shitf(pi^e(U_d,0^e), z^e)
            pos_pi = 0;
            for (size_t k=0; k<word_K-1; k++) {
                shift_u[k] = rotl16((uint16_t) extract_sixteen_bits(U[e][d][pi[pos_pi]/4], pi[pos_pi]%4), extract_four_bits(z[e][pos_pi/16], pos_pi%16));
                shift_u[k] ^= (u64) rotl16((uint16_t) extract_sixteen_bits(U[e][d][pi[(pos_pi+1)]/4], pi[(pos_pi+1)]%4), extract_four_bits(z[e][(pos_pi+1)/16], (pos_pi+1)%16)) << 16;
                shift_u[k] ^= (u64) rotl16((uint16_t) extract_sixteen_bits(U[e][d][pi[(pos_pi+2)]/4], pi[(pos_pi+2)]%4), extract_four_bits(z[e][(pos_pi+2)/16], (pos_pi+2)%16)) << 32;
                shift_u[k] ^= (u64) rotl16((uint16_t) extract_sixteen_bits(U[e][d][pi[(pos_pi+3)]/4], pi[(pos_pi+3)]%4), extract_four_bits(z[e][(pos_pi+3)/16], (pos_pi+3)%16)) << 48;

                pos_pi += 4;
            }
            shift_u[word_K-1] = rotl16((uint16_t) extract_sixteen_bits(U[e][d][pi[pos_pi]/4], pi[pos_pi]%4), 
                            extract_four_bits(z[e][pos_pi/16], pos_pi%16));
            
            mul_mat_vec(y_d[0], H, shift_u);
            // y_d1^e <- y_d0^e xor y
            for (size_t k=0; k<word_k; k++) {
                y_d[1][k] = y_d[0][k] ^ y[k];
            }
            //
            pos_pi = 0;
            for (size_t k=0; k<word_x-1; k++) {
                for (size_t pos=0; pos<BS; pos++) { 
                    // z_d0^e <- X_d^e - pi^e(R_d0^e) mod bs
                    z_d[0][k] = (z_d[0][k] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(X[e][d][k], pos) 
                            - extract_four_bits(R[e][d][pi[pos_pi]/16], pi[pos_pi]%16)) % BS) << (4*pos));
                    pos_pi++;
                }
            }
            for (size_t pos=0; pos<13; pos++) { 
                z_d[0][word_x-1] = (z_d[0][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(X[e][d][word_x-1], pos) 
                        - extract_four_bits(R[e][d][pi[pos_pi]/16], pi[pos_pi]%16)) % BS) << (4*pos));
                pos_pi++;
            }
            for (size_t k=0; k<word_x-1; k++) {
                for (size_t pos=0; pos<BS; pos++) { 
                    // z_d1^e <- z^e - z_d0^e  mod bs
                    z_d[1][k] = (z_d[1][k] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(z[e][k], pos) 
                            - extract_four_bits(z_d[0][k], pos)) % BS) << (4*pos));
                }
            }
            for (size_t pos=0; pos<13; pos++) { 
                z_d[1][word_x-1] = (z_d[1][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(z[e][word_x-1], pos) 
                            - extract_four_bits(z_d[0][word_x-1], pos)) % BS) << (4*pos));
            }
            
            EVP_DigestUpdate(h2, y_d[0], word_k*8);  
            EVP_DigestUpdate(h2, y_d[1], word_k*8);
            EVP_DigestUpdate(h2, z_d[0], word_x*8);
            EVP_DigestUpdate(h2, z_d[1], word_x*8);
        }
    }

    //// Phase 4
    // h2 <- H2(m, salt, h1, (y^e_db, z^e_db))
    // Set (b_1^e, ..., b_D^e) <- PRG_2(h2)
    // Let i^e <- sum_{d=1}^D b_d^e 2^{d-1}
    EVP_DigestFinal_ex(h2, buff_h2, &resultLen);
    memcpy(sigma->h2, buff_h2, SHA256_DIGEST_LENGTH);
    
    u64 out_copath[tau][D][2];
    size_t i_e;
    size_t lenb = 0;
    for (size_t e=0; e<tau; e++) {
        i_e = 0;
        for (size_t d=0; d<D; d++) {
            if (buff_h2[lenb/8] & 1) {
                i_e ^= 1<<d;
            }
            buff_h2[lenb/8] >>= 1;
            lenb += 1;
        }
        #ifdef PPRF_SHA
        compute_copath_sha(out_copath[e], seed_e[e], i_e, e, K_0, K_1);
        #else
        compute_copath(out_copath[e], seed_e[e], i_e, key_schedule_k0, key_schedule_k1);
        #endif
        

        if (i_e != (n_par-1)) {
            memcpy(&sigma->aux.x[e], &x[e][0], word_x*8);
            memcpy(&sigma->aux.u_n[e], &u[e][0], word_K*8);
            memset(&sigma->com_n[e], 0, SHA256_DIGEST_LENGTH);
            memcpy(&sigma->com[e], &log_com[e][i_e][0], lambda/8);
        }
    }

    //// Phase 5
    // Output sigma = (salt, h1, h2, z^e, CoPath_salt(i^e, seed^e), com_i^e^e, aux_n^e)
    memcpy(&sigma->sig_salt.K_0, K_0, lambda/8);
    memcpy(&sigma->sig_salt.K_1, K_1, lambda/8);
    memcpy(&sigma->copath, out_copath, tau*D*8*2);
    memcpy(&sigma->z, z, tau*word_x*8);
    
    return 0;
}

// Signature algorithm of the signature scheme 
int folding_buckets (__m128i key_schedule_k0[11]) {
    //// Phase 1
    u64 X[tau][D][word_x], R[tau][D][word_x], U[tau][D][word_K];
    u64 x[tau][word_x], u[tau][word_K], r[tau][word_x];
    u64 r_n[tau][word_x];
    u64 ue[word_K];

    for (size_t e=0; e<tau; e++) {
        for (size_t d=0; d<D; d++) {
            for (size_t k=0; k<word_x; k++) {
                X[e][d][k] = 0;
                R[e][d][k] = 0;
            }
            for (size_t k=0; k<word_K; k++) {
                U[e][d][k] = 0;      
            }
        }
        for (size_t k=0; k<word_x; k++) {
            x[e][k] = 0;
            r[e][k] = 0;
            r_n[e][k] = 0;
        }
        for (size_t k=0; k<word_K; k++) {
            u[e][k] = 0;
        }
    }
    
    u64 sampled_data[((lambda + big_K + 2 * log_bs * w)/64)+ 8];
    u64 extracted_pos = 0;


    u64 seed_i[word_lambda] = {9322567240304977384ULL, 12730924059968730833ULL};
    u64 scal_a = 16664525ULL;
    u64 scal_b = 1013904223ULL;
    u64 modulus = (u64)-1;
    for (size_t e=0; e<tau; e++) {
        for (size_t i=0; i<n_par-1; i++) {
            seed_i[0] = (scal_a * seed_i[0] + scal_b) % modulus;
            seed_i[1] = (scal_a * seed_i[1] + scal_b) % modulus;
            // (x_i^e, r_i^e, u_i^e, com_i^e) <- PRG(seed_i^e)
            sample_random_from_seed(seed_i, sampled_data, lambda + big_K + 2*(log_bs * w)+ 4*64, key_schedule_k0);
            // x_n^e <- x_n^e - x_i^e mod bs, u_n^e <- u_n^e xor u_i^e, and r^e <- r^e+r_i^e mod bs
            extracted_pos = 0;
            for (size_t k=0; k<word_x-1; k++) {
                for (size_t pos=0; pos<BS; pos++) {
                    x[e][k] = (x[e][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(x[e][k], pos) 
                                - extract_four_bits(sampled_data[k], pos)) % BS) << (4*pos));
                }
            }
            for (size_t pos=0; pos<13; pos++) {
                x[e][word_x-1] = (x[e][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(x[e][word_x-1], pos) 
                            - extract_four_bits(sampled_data[word_x-1], pos)) % BS) << (4*pos));
            }
            for (size_t k=0; k<word_x-1; k++) {
                for (size_t pos=0; pos<BS; pos++) {
                    r[e][k] = (r[e][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(r[e][k], pos) 
                                + extract_four_bits(sampled_data[word_x+k], pos)) % BS) << (4*pos));
                }
            }
            for (size_t pos=0; pos<13; pos++) {
                r[e][word_x-1] = (r[e][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(r[e][word_x-1], pos) 
                            + extract_four_bits(sampled_data[word_x+word_x-1], pos)) % BS) << (4*pos));
            }
            extracted_pos += 2*word_x;
            for (size_t k=0; k<word_K; k++) {
                u[e][k] ^= sampled_data[extracted_pos + k];
            }
            u[e][word_K-1] = u[e][word_K-1] & 0xFFFF;
            extracted_pos += word_K;
            
            // For all d <= D such that i[d] = 0, set
            for (size_t d=0; d<D; d++) {
                if (((i>>d) & 1) == 0) {
                    // X_d,0^e <- X_d,0^e + x_i^e mod bs
                    // R_d,0^e <- R_d,0^e + r_i^e mod bs
                    for (size_t k=0; k<word_x-1; k++) {
                        for (size_t pos=0; pos<BS; pos++) {
                            X[e][d][k] = (X[e][d][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[k], pos) 
                                        + extract_four_bits(X[e][d][k], pos)) % BS) << (4*pos));                           
                        }
                    }
                    for (size_t pos=0; pos<13; pos++) {
                            X[e][d][word_x-1] = (X[e][d][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x-1], pos) 
                                        + extract_four_bits(X[e][d][word_x-1], pos)) % BS) << (4*pos));
                    }
                    for (size_t k=0; k<word_x-1; k++) {
                        for (size_t pos=0; pos<BS; pos++) {
                            R[e][d][k] = (R[e][d][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x + k], pos) 
                                        + extract_four_bits(R[e][d][k], pos)) % BS) << (4*pos));
                           
                        }
                    }
                    for (size_t pos=0; pos<13; pos++) {
                            R[e][d][word_x-1] = (R[e][d][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x+word_x-1], pos) 
                                        + extract_four_bits(R[e][d][word_x-1], pos)) % BS) << (4*pos));
                    }
                    // U_d,0^e <- U_d,0^e xor u_i^e
                    for (size_t k=0; k<word_K; k++) {
                        U[e][d][k] ^= sampled_data[2*word_x + k];
                    }
                    U[e][d][word_K-1] = U[e][d][word_K-1] & 0xFFFF;
                }
            }
        }
        // Compute r_n^e <- PRG(seed_n^e)
        seed_i[0] = 11906471602029514385ULL;
        seed_i[1] = 13681627746222550614ULL;
        sample_random_from_seed(seed_i, r_n[e], log_bs*w, key_schedule_k0);
        // r^e <- r^e + r_n^e mod bs, u^e <- Expand(r^e), and u_n^e<-u_n^e xor u^e
        for (size_t k=0; k<word_x-1; k++) {
            for (size_t pos=0; pos<BS; pos++) {
                r[e][k] = (r[e][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(r[e][k], pos) 
                            + extract_four_bits(r_n[e][k], pos)) % BS) << (4*pos));
            }
        }
        for (size_t pos=0; pos<13; pos++) {
            r[e][word_x-1] = (r[e][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(r[e][word_x-1], pos) 
                        + extract_four_bits(r_n[e][word_x-1], pos)) % BS) << (4*pos));
        }
        for (size_t k=0; k<word_K; k++) {
            ue[k] = 0;
        }
        expand(r[e], ue);
        for (size_t k=0; k<word_K; k++) {
            u[e][k] ^= ue[k];
        }
    }
   return 0;
}


int crypto_sign_open (u64 *m, 
           size_t m_len,
           u64 H[word_H], 
           u64 y[word_k], 
           __m128i key_schedule_k0[11], 
           __m128i key_schedule_k1[11],
           sig_t *signature 
) {

     //// Initialization
    unsigned int resultLen = 0;
    unsigned char buff_h1[SHA256_DIGEST_LENGTH];
    unsigned char buff_h2[SHA256_DIGEST_LENGTH];
    EVP_MD_CTX* h1 = EVP_MD_CTX_new();
    EVP_DigestInit_ex(h1, EVP_sha256(), NULL);
    EVP_MD_CTX* h2 = EVP_MD_CTX_new();
    EVP_DigestInit_ex(h2, EVP_sha256(), NULL);

    EVP_DigestUpdate(h1, m, m_len);
    EVP_DigestUpdate(h1, signature->sig_salt.K_0, lambda/8);
    EVP_DigestUpdate(h1, signature->sig_salt.K_1, lambda/8);
    
    EVP_DigestUpdate(h2, m, m_len);
    EVP_DigestUpdate(h2, signature->sig_salt.K_0, lambda/8);
    EVP_DigestUpdate(h2, signature->sig_salt.K_1, lambda/8);

 
    // Recompute (b1^e, ..., b_D^e) using h2 abd define the i^e
    size_t i_e;
    size_t lenb = 0;
    unsigned char tmp_sig_h2[SHA256_DIGEST_LENGTH];

    memcpy(tmp_sig_h2, signature->h2, SHA256_DIGEST_LENGTH);

    u64 seed_i[word_lambda];
    u64 X[2][tau][D][word_x];
    u64 R[2][tau][D][word_x];
    u64 U[2][tau][D][word_K]; 

    for (size_t e=0; e<tau; e++) {
        for (size_t d=0; d<D; d++) {
            for (size_t k=0; k<word_x; k++) {
                X[0][e][d][k] = 0;
                X[1][e][d][k] = 0;
                R[0][e][d][k] = 0;
                R[1][e][d][k] = 0;
            }
            for (size_t k=0; k<word_K; k++) {
                U[0][e][d][k] = 0;
                U[1][e][d][k] = 0;
            }
        }
    }

    u64 sampled_data[((lambda + big_K + 2 * log_bs * w)/64)+ 8];
    size_t tab_ie[tau];
    u64 r_n[tau][word_x];
    for (size_t e=0; e<tau; e++) {
        for (size_t k=0; k<word_x; k++) {
            r_n[e][k] = 0;
        }
    }

    for (size_t e=0; e<tau; e++) {
        // recompute i_e
        i_e = 0;
        for (size_t d=0; d<D; d++) {
            if (tmp_sig_h2[lenb/8] & 1) {
                i_e ^= 1<<d;
            }
            tmp_sig_h2[lenb/8] >>= 1;
            lenb += 1;
        }
        tab_ie[e] = i_e;
        for (size_t i=0; i<n_par; i++) {
            if (i == i_e) {
                if (i == (n_par-1)) {
                    EVP_DigestUpdate(h1, signature->com_n[e], SHA256_DIGEST_LENGTH);
                }
                else {
                    EVP_DigestUpdate(h1, signature->com[e], lambda/8);
                }
                continue;
            }
            // recompute seed_i_e^e
            #ifdef PPRF_SHA
            recover_seed_from_copath_sha(seed_i, signature->copath[e], i_e, i, e, signature->sig_salt.K_0, signature->sig_salt.K_1);
            #else
            recover_seed_from_copath(seed_i, signature->copath[e], i_e, i, key_schedule_k0, key_schedule_k1);
            #endif
            sample_random_from_seed(seed_i, sampled_data, lambda + big_K + 2*(log_bs * w) + 4*64, key_schedule_k0);
            // For all d <= D such that i[d] = b_i, set
            if (i != n_par-1) {
                for (size_t d=0; d<D; d++) {
                    size_t b_e = 1 - ((i_e>>d) & 1);
                    if (((i>>d) & 1) != ((i_e>>d) & 1)) {
                        // X_d,b^e <- X_d,b^e + x_i^e mod bs
                        // R_d,b^e <- R_d,b^e + r_i^e mod bs
                        for (size_t k=0; k<word_x-1; k++) {
                            for (size_t pos=0; pos<BS; pos++) {
                                X[b_e][e][d][k] = (X[b_e][e][d][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[k], pos) 
                                        + extract_four_bits(X[b_e][e][d][k], pos)) % BS) << (4*pos));
                            }
                        }
                        for (size_t pos=0; pos<13; pos++) {
                                X[b_e][e][d][word_x-1] = (X[b_e][e][d][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x-1], pos) 
                                        + extract_four_bits(X[b_e][e][d][word_x-1], pos)) % BS) << (4*pos));
                        }
                        for (size_t k=0; k<word_x-1; k++) {
                            for (size_t pos=0; pos<BS; pos++) {
                                R[b_e][e][d][k] = (R[b_e][e][d][k] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x + k], pos) 
                                        + extract_four_bits(R[b_e][e][d][k], pos)) % BS) << (4*pos));
                            }
                        }
                        for (size_t pos=0; pos<13; pos++) {
                                R[b_e][e][d][word_x-1] = (R[b_e][e][d][word_x-1] & (~(((u64)(0xF))<<(pos*4))))  ^ (u64)(((extract_four_bits(sampled_data[word_x + word_x-1], pos) 
                                        + extract_four_bits(R[b_e][e][d][word_x-1], pos)) % BS) << (4*pos));
                               
                        }
                        // U_d,0^e <- U_d,0^e xor u_i^e
                        for (size_t k=0; k<word_K; k++) {
                            U[b_e][e][d][k] ^= sampled_data[2*word_x + k]; 
                        }
                        U[b_e][e][d][word_K-1] = U[b_e][e][d][word_K-1] & 0xFFFF;
                        
                    } 
                }
                 // h(com_i^e)
                EVP_DigestUpdate(h1, &sampled_data[2*word_x+word_K], lambda/8);
            }
            else if (i==n_par-1) {
                // Compute r_n^e <- PRG(seed_n^e)
                sample_random_from_seed(seed_i, r_n[e], (size_t) log_bs*w, key_schedule_k0);
                for (size_t d=0; d<D; d++) {
                    size_t b_e = 1 - ((i_e>>d) & 1);
                    if (((i>>d) & 1) != ((i_e>>d) & 1)) {
                        // X_d,b^e <- X_d,b^e + x_i^e mod bs
                        // R_d,b^e <- R_d,b^e + r_i^e mod bs
                        for (size_t k=0; k<word_x-1; k++) {
                            for (size_t pos=0; pos<BS; pos++) {
                                X[b_e][e][d][k] = (X[b_e][e][d][k] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(signature->aux.x[e][k], pos) 
                                            + extract_four_bits(X[b_e][e][d][k], pos)) % BS) << (4*pos));
                            }
                        }
                        for (size_t pos=0; pos<13; pos++) {
                                X[b_e][e][d][word_x-1] = (X[b_e][e][d][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(signature->aux.x[e][word_x-1], pos) 
                                            + extract_four_bits(X[b_e][e][d][word_x-1], pos)) % BS) << (4*pos));
                        }
                        for (size_t k=0; k<word_x-1; k++) {
                            for (size_t pos=0; pos<BS; pos++) {
                                R[b_e][e][d][k] = (R[b_e][e][d][k] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(r_n[e][k], pos) 
                                            + extract_four_bits(R[b_e][e][d][k], pos)) % BS) << (4*pos));
                            }
                        }
                        for (size_t pos=0; pos<13; pos++) {
                                R[b_e][e][d][word_x-1] = (R[b_e][e][d][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(r_n[e][word_x-1], pos) 
                                            + extract_four_bits(R[b_e][e][d][word_x-1], pos)) % BS) << (4*pos));
                        }
                        // U_d,0^e <- U_d,0^e xor u_i^e
                        for (size_t k=0; k<word_K; k++) {
                            U[b_e][e][d][k] ^= signature->aux.u_n[e][k];
                        }
                    } 
                }
                // com_n^e <-H(state_n^e)
                unsigned char buff_h[SHA256_DIGEST_LENGTH];
                EVP_MD_CTX* h = EVP_MD_CTX_new();
                EVP_DigestInit_ex(h, EVP_sha256(), NULL);
                EVP_DigestUpdate(h, &signature->aux.x[e][0], word_x*8);
                EVP_DigestUpdate(h, &signature->aux.u_n[e][0], word_K*8);
                EVP_DigestUpdate(h, seed_i, lambda/8);
                EVP_DigestFinal_ex(h, buff_h, &resultLen);
                // Update h1 with com_n^e
                EVP_DigestUpdate(h1, buff_h, SHA256_DIGEST_LENGTH);
            }
        }
        
    }

    EVP_DigestFinal_ex(h1, buff_h1, &resultLen);
    int res_h1 = memcmp(buff_h1, signature->h1, SHA256_DIGEST_LENGTH);
    EVP_DigestUpdate(h2, buff_h1, SHA256_DIGEST_LENGTH);

    // Recompute pi^e using h1
    uint8_t pi[w];
    uint8_t rand_pi[tau][w];
    sample_random_from_seed((uint64_t*)&buff_h1[0], (uint64_t*) &rand_pi[0][0], tau*w*8, key_schedule_k0);

    //// Phase 3
    // For each iteration e in tau
    u64 z_d[2][word_x];
    u64 y_d[2][word_k];
    u64 shift_u[word_K];
    size_t pos_pi;

    
    size_t b_e = 0;
    for (size_t e=0; e<tau; e++) {
        // generate random permutation using Fisher-Yates
        for (size_t i=0; i<w; i++) {
            pi[i] = i;
        }
        for (size_t i=w-1; i<w; --i) {
            size_t j = rand_pi[e][i] % (i+1);
            size_t tmp = pi[i];
            pi[i] = pi[j];
            pi[j] = tmp;
        }

        for (size_t d=0; d<D; d++) {
            b_e = 1 - ((tab_ie[e] >> d) & 1);
            // y_d0^e <- H * Shitf(pi^e(U_d,0^e), z^e)
            pos_pi = 0;
            for (size_t k=0; k<word_K-1; k++) {
                shift_u[k] = rotl16((uint16_t) extract_sixteen_bits(U[b_e][e][d][pi[pos_pi]/4], pi[pos_pi]%4), extract_four_bits(signature->z[e][pos_pi/16], pos_pi%16));
                shift_u[k] ^= (u64) rotl16((uint16_t) extract_sixteen_bits(U[b_e][e][d][pi[(pos_pi+1)]/4], pi[(pos_pi+1)]%4), extract_four_bits(signature->z[e][(pos_pi+1)/16], (pos_pi+1)%16)) << 16;
                shift_u[k] ^= (u64) rotl16((uint16_t) extract_sixteen_bits(U[b_e][e][d][pi[(pos_pi+2)]/4], pi[(pos_pi+2)]%4), extract_four_bits(signature->z[e][(pos_pi+2)/16], (pos_pi+2)%16)) << 32;
                shift_u[k] ^= (u64) rotl16((uint16_t) extract_sixteen_bits(U[b_e][e][d][pi[(pos_pi+3)]/4], pi[(pos_pi+3)]%4), extract_four_bits(signature->z[e][(pos_pi+3)/16], (pos_pi+3)%16)) << 48;
                pos_pi += 4;   
            }
            shift_u[word_K-1] = rotl16((uint16_t)extract_sixteen_bits(U[b_e][e][d][pi[pos_pi]/4], pi[pos_pi]%4), extract_four_bits(signature->z[e][pos_pi/16], pos_pi%16));

            mul_mat_vec(y_d[b_e], H, shift_u);
            // y_d1^e <- y_d0^e xor y
            for (size_t k=0; k<word_k; k++) {
                y_d[1-b_e][k] = y_d[b_e][k] ^ y[k]; 
            }
            //
            pos_pi = 0;
            for (size_t k=0; k<word_x-1; k++) {
                for (size_t pos=0; pos<BS; pos++) { 
                    // z_dbe^e <- X_d^ebe - pi^e(R_dbe^e) mod bs
                    z_d[b_e][k] = (z_d[b_e][k] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(X[b_e][e][d][k], pos) 
                            - extract_four_bits(R[b_e][e][d][pi[pos_pi]/16], pi[pos_pi]%16)) % BS) << (4*pos));
                    pos_pi++;
                }
            }
            for (size_t pos=0; pos<13; pos++) { 
                z_d[b_e][word_x-1] = (z_d[b_e][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(X[b_e][e][d][word_x-1], pos) 
                        - extract_four_bits(R[b_e][e][d][pi[pos_pi]/16], pi[pos_pi]%16)) % BS) << (4*pos));
                pos_pi++;
            }
            for (size_t k=0; k<word_x-1; k++) {
                for (size_t pos=0; pos<BS; pos++) { 
                    // z_d1^e <- z^e - z_dbe^e  mod bs
                    z_d[1-b_e][k] = (z_d[1-b_e][k] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(signature->z[e][k], pos) 
                            - extract_four_bits(z_d[b_e][k], pos)) % BS) << (4*pos));
                }
            }
            for (size_t pos=0; pos<13; pos++) { 
                z_d[1-b_e][word_x-1] = (z_d[1-b_e][word_x-1] & (~(((u64)(0xF))<<(pos*4)))) ^ (u64)(((extract_four_bits(signature->z[e][word_x-1], pos) 
                        - extract_four_bits(z_d[b_e][word_x-1], pos)) % BS) << (4*pos));
            }
            EVP_DigestUpdate(h2, y_d[0], word_k*8);
            EVP_DigestUpdate(h2, y_d[1], word_k*8);
            EVP_DigestUpdate(h2, z_d[0], word_x*8);
            EVP_DigestUpdate(h2, z_d[1], word_x*8);
        }
    }

    EVP_DigestFinal_ex(h2, buff_h2, &resultLen);
    int res_h2 = memcmp(buff_h2, signature->h2, SHA256_DIGEST_LENGTH);
    
    #ifdef PPRF_SHA
    return 0;
    #else
    if (res_h2 != 0) {
        printf("sign:    ");
        for (size_t v=0; v<SHA256_DIGEST_LENGTH; v++) {
            printf("%X", signature->h2[v]);
        }
        printf("\nverify:  ");
        for (size_t v=0; v<SHA256_DIGEST_LENGTH; v++) {
            printf("%X", buff_h2[v]);
        }
        printf("\n");
    }
    return (res_h1 | res_h2);
    #endif
}
