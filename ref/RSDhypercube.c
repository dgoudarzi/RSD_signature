/**
 * @file
 * This file contains the RSD signatures implementation.
 *     Copyright (C) 2024  Antoine Joux, Dahmun Goudarzi 
 *
 * D  |   tau
 * ---|------
 * 8  |   16
 * 9  |   15
 * 10 |   13
 * 11 |   12
 * 12 |   11
 * 13 |   10
 * 15 |   9
 * 16 |   8
 *
 * To compile, test and time all versions, simply run maketest.sh
 **/
#include "RSDinclude.h"

// Choose DBLing function
#ifdef CLASSICAL_SHA_TREE
#ifdef NEED_XOR_PRESERVE
#define AES_DBL sha256_dbl_enc_delta_preserve
#define AES_EXPAND sha256_expand_1cpy
#else
#define AES_DBL sha256_dbl_enc
#define AES_EXPAND sha256_expand
#endif
#else
#ifdef NEED_XOR_PRESERVE
#define AES_DBL aes128_dbl_enc_delta_preserve
#define AES_EXPAND aes128_expand_1cpy
#else
#define AES_DBL aes128_dbl_enc
#define AES_EXPAND aes128_expand
#endif
#endif

/*! \brief Multiply Matrix with Sparse Vector
 */
void SparseMatMul(u64 *Synd, u64 *SecX, u64 *PubMat)
{
  int i, j, pos;
  
  for (j = 0; j < 2 * CodeRepinU128; j++) 
    Synd[j] = 0;

  for (i = 0; i < NumWindows; i++)
  {
    pos = (SecX[i / 16] >> (4 * (i % 16))) & 15; /* In-window position on 4 bits to accommodate 8 and 16 */
    for (j = 0; j < 2 * CodeRepinU128; j++)
      Synd[j] ^= PubMat[2 * (i * WindowModulus + pos) * CodeRepinU128 + j];
  }
}

/*! \brief Multiply Matrix with Dense Vector
 */
void DenseMatMul(u64 *Synd, u64 *VecX, u64 *PubMat)
{
  int i, j, pos;
  u64 mask;

  mask = -(VecX[0] & 1);

#pragma GCC unroll 4
#pragma GCC ivdep
  for (j = 0; j < 2 * CodeRepinU128; j++)
    Synd[j] = PubMat[2 * i * CodeRepinU128 + j] & mask;

  for (i = 1; i < CodeLength; i++)
  {
    mask = -((VecX[i / 64] >> (i % 64)) & 1);
#pragma GCC unroll 4
#pragma GCC ivdep
    for (j = 0; j < 2 * CodeRepinU128; j++)
      Synd[j] ^= PubMat[2 * i * CodeRepinU128 + j] & mask;
  }
}

/*! \brief Multiply Matrix with Sparse Vector adding result to Accumulator
 */
void DenseMatMulAdd(u64 *Synd, u64 *VecX, u64 *PubMat)
{
  int i, j, pos;

  for (i = 0; i < CodeLength; i++)
  {
    for (i = 0; i < CodeLength; i++)
    {
      u64 mask;
      mask = -((VecX[i / 64] >> (i % 64)) & 1);
#pragma GCC unroll 4
#pragma GCC ivdep
      for (j = 0; j < 2 * CodeRepinU128; j++)
        Synd[j] ^= PubMat[2 * i * CodeRepinU128 + j] & mask;
    }
  }
}

/*! \brief Expand Public Matrix from Seed
 */
void GeneratePubMat(u64 *PubMat, u64 *PubSeed)
{
  aessubkeytype key_schedule[11];
  int i, j, pos, zoffset;
  u64 zmask;

  aes128_expand_key(key_schedule, PubSeed);

  for (i = 0; i < CodeLength; i++)
  {
    for (j = 0; j < 2 * CodeRepinU128; j++)
    {
      pos = 2 * i * CodeRepinU128 + j;
      PubMat[pos] = pos;
    }
  }
  aes128_enc(PubMat, PubMat, CodeLength * CodeRepinU128, key_schedule);

  zoffset = CodeNumChecks % 128;
  if (zoffset != 0)
  {
    if (zoffset < 64)
    {
      zmask = (1UL << zoffset) - 1;
      for (i = 0; i < CodeLength; i++)
      {
        pos = 2 * i * CodeRepinU128 + 2 * CodeRepinU128 - 2;
        PubMat[pos] &= zmask;
        PubMat[pos + 1] = 0;
      }
    }
    else
    {
      zmask = (1UL << (zoffset - 64)) - 1;
      for (i = 0; i < CodeLength; i++)
      {
        pos = 2 * i * CodeRepinU128 + 2 * CodeRepinU128 - 2;
        PubMat[pos + 1] &= zmask;
      }
    }
  }
}

/*! \brief Create public key from prescribed secret key
 */
void CreateKeyFromSeed(u64 Seed[4], u64 *SecX, u64 *Synd, u64 *PubMat, u64 *PubSeed)
{
  SHA256_CTX hash_ctx;
  u64 inner_hash[4];
  u64 Tmp[5];
  u64 Tmp2[5];
  aessubkeytype key_schedule[11];
  int i, j, pos;

  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, (uchar *)Seed, 32);
  SHA256Final((uchar *)inner_hash, &hash_ctx);

  aes128_expand_key(key_schedule, inner_hash);

  for (i = 0; i < 2 * SecretRepinU128; i++)
  {
    SecX[i] = i;
  }
  aes128_enc(SecX, SecX, SecretRepinU128, key_schedule);
  for (i = 0; i < 2 * (SecretRepinU128 - 1); i++)
  {
    SecX[i] &= SecretInnerMask;
  }
  SecX[i] &= SecretInnerMask;    // Specific to shape of 217 windows mod 8
  SecX[i + 1] &= SecretLastMask; // Specific to shape of 217 windows mod 8

  PubSeed[0] = inner_hash[2];
  PubSeed[1] = inner_hash[3];

  GeneratePubMat(PubMat, PubSeed);
  SparseMatMul(Synd, SecX, PubMat);
}

/*! \brief Perform the key schedule for all AES keys used in the scheme
 */
void FillAESKey(u64 *AES_keys, aessubkeytype *key_schedules, int nbkeys)
{
  int i;
  aessubkeytype local_key_schedule[11];

  for (i = 0; i < 2; i++)
    AES_keys[i + 2] = 0;

  // Encrypt K0 with cstant key=0 to get K1, then K2,...
  aes128_expand_key(local_key_schedule, AES_keys + 2);
  aes128_enc(AES_keys, AES_keys + 2, nbkeys - 1, local_key_schedule);

  for (i = 0; i < nbkeys; i++)
  {
    aes128_expand_key(key_schedules + 11 * i, AES_keys + 2 * i);
  }
}

/*! \brief Initialize keys and pseudo-random values from Message and Global seed
 */
void InitKeyMaterial(uchar *messhash, uchar *randomness,
                     u64 *AES_keys, aessubkeytype *key_schedules, int nbkeys,
                     uchar *inner_hash, int nbrand)
{
  SHA256_CTX hash_ctx;
  uchar local_hash[32];
  int i;

  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, messhash, 32);
  SHA256Update(&hash_ctx, randomness, 32);
  SHA256Final(local_hash, &hash_ctx);

  for (i = 0; i < 2; i++)
    AES_keys[i] = 0;
  for (i = 0; i < 16; i += 2)
  {

    AES_keys[0] = (AES_keys[0] << 8) | local_hash[i];
    AES_keys[1] = (AES_keys[1] << 8) | local_hash[i + 1];
  }

  FillAESKey(AES_keys, key_schedules, nbkeys);

  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, local_hash, 32);
  SHA256Final(inner_hash, &hash_ctx);
  for (i = 0; i < nbrand - 1; i++)
  {
    SHA256Init(&hash_ctx);
    SHA256Update(&hash_ctx, inner_hash + 32 * i, 32);
    SHA256Final(inner_hash + 32 * (i + 1), &hash_ctx);
  }
}

/*! \brief Prepare the first level nodes of all PPRF trees
 */
void InitTrees(int NumRounds, u64 *FixedDifference, int needfixed,
               uchar *randomness, u64 *FullTrees)
{
  int i, r;
  u64 TreeOffset;

  TreeOffset = 4 * (LEAVES_SIZE - 2 * NumRounds);
  for (r = 0; r < 4 * NumRounds; r++)
    FullTrees[TreeOffset + r] = 0;
  for (r = 0; r < NumRounds; r++)
  {
    TreeOffset = 4 * (LEAVES_SIZE - 2 * NumRounds) + 4 * r;
    for (i = 0; i < 16; i += 2)
    {
      FullTrees[TreeOffset] = (FullTrees[TreeOffset] << 8) | (randomness[i + 32 * r]);
      FullTrees[TreeOffset + 1] = (FullTrees[TreeOffset + 1] << 8) | (randomness[i + 32 * r + 1]);
    }
    if (needfixed)
    {
      FullTrees[TreeOffset + 2] = FullTrees[TreeOffset] ^ FixedDifference[0];
      FullTrees[TreeOffset + 3] = FullTrees[TreeOffset + 1] ^ FixedDifference[1];
    }
    else
    {
      for (i = 0; i < 16; i += 2)
      {
        FullTrees[TreeOffset + 2] = (FullTrees[TreeOffset + 2] << 8) | (randomness[i + 32 * r + 16]);
        FullTrees[TreeOffset + 3] = (FullTrees[TreeOffset + 3] << 8) | (randomness[i + 32 * r + 17]);
      }
    }
  }
}

/*! \brief Fill all PPRF trees and perform expansion of the leaves
 */
void GenerateTrees(int NumRounds, int lvl_offset,
                   aessubkeytype *key_schedules,
                   u64 *FullTrees, u64 *ExpandedTree, int *Queries)
{
  u64 CurrentSize;
  int lvl, i, r, pos;

  // Tree
  for (CurrentSize = 2 * NumRounds, lvl = 1; CurrentSize < LEAVES_SIZE; CurrentSize <<= 1, lvl++)
  {
    u64 SaveVals[2 * NbRounds];
    u64 *SaveTree, SavePos;
    int r;
    if (Queries != NULL)
    {
      SaveTree = FullTrees + 4 * (LEAVES_SIZE - 2 * CurrentSize);
      for (r = 0; r < NbRounds; r++)
      {
        SavePos = ((Queries[r] >> (LogPlayers - 1 - lvl)) ^ 1) + r * (2 * CurrentSize / NbRounds);
        SaveVals[2 * r] = SaveTree[2 * SavePos];
        SaveVals[2 * r + 1] = SaveTree[2 * SavePos + 1];
      }
    }
    AES_DBL(FullTrees + 4 * (LEAVES_SIZE - CurrentSize),
            FullTrees + 4 * (LEAVES_SIZE - 2 * CurrentSize),
            CurrentSize, key_schedules + 2 * (lvl + lvl_offset) * 11,
            key_schedules + (2 * (lvl + lvl_offset) + 1) * 11);
    if (Queries != NULL)
    {
      for (r = 0; r < NbRounds; r++)
      {
        SavePos = ((Queries[r] >> (LogPlayers - 1 - lvl)) ^ 1) + r * (2 * CurrentSize / NbRounds);
        SaveTree[2 * SavePos] = SaveVals[2 * r];
        SaveTree[2 * SavePos + 1] = SaveVals[2 * r + 1];
      }
    }
  }

  // Leaves expansion
  AES_EXPAND(FullTrees, ExpandedTree,
             EXPAND_RATIO, LEAVES_SIZE,
             key_schedules + (2 * (lvl + lvl_offset) + 1) * 11);

  // Remove high bits of modular window part -- Needed for Folding
  for (i = 0; i < LEAVES_SIZE; i++)
  {
    for (pos = 0; pos < EXPAND_RATIO; pos++)
      ExpandedTree[2 * i * EXPAND_RATIO + pos] &= SecretInnerMask;
  }
}

void ClearUnused(int NumRounds, u64 *FoldedExpTree, u64 *Offsets)
{
  int r, i, pos;

  for (r = 0; r < NumRounds; r++)
  {
    Offsets[2 * EXPAND_RATIO * r + 13] &= SecretLastMask;
    Offsets[2 * EXPAND_RATIO * r + 14 + 13] &= SecretLastMask;
    Offsets[2 * EXPAND_RATIO * r + 28 + 27] &= 0xff;

    for (i = 0; i < 2 * LogPlayers; i++)
    {
      FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + 13] &= SecretLastMask;
      FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + 14 + 13] &= SecretLastMask;
      FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + 28 + 27] &= 0xff;
#ifdef COMPRESSION
      for (pos = EXPAND_RATIO; pos < 2 * EXPAND_RATIO; pos++)
        FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + pos] &= CompressMask;
#endif
    }
  }
}

#ifdef CLASSICAL_FOLDING
void FoldTrees(int NumRounds, u64 *ExpandedTree, u64 *FoldedExpTree, u64 *Offsets)
#else
void SlowFoldTrees(int NumRounds, u64 *ExpandedTree, u64 *FoldedExpTree, u64 *Offsets)
#endif
{
  int i, pos, r;
  for (i = 0; i < 4 * EXPAND_RATIO * (NumRounds * LogPlayers); i++)
    FoldedExpTree[i] = 0;

  for (r = 0; r < NumRounds; r++)
  {
    for (i = 0; i < NbPlayers; i++)
    {
      for (int bitpos = 0; bitpos < LogPlayers - 1; bitpos++)
      {
        if (((i >> bitpos) & 1) == 0)
        {
          // Half modular, half binary
#pragma GCC unroll 4
#pragma GCC ivdep
          for (pos = 0; pos < EXPAND_RATIO; pos++)
          {
            FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos] +=
                ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos];
            FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos] &=
                SecretInnerMask;
            // Get rid of carries i.e. reduces mod 8 [Does not work for 16]
          }
          for (; pos < 2 * EXPAND_RATIO; pos++)
          {
            FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos] ^=
                ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos];
          }
        }
      }
    }
  }

  for (r = 0; r < NumRounds; r++)
  {
    for (i = 0; i < NbPlayers / 2; i++)
    {
#pragma GCC unroll 4
#pragma GCC ivdep
      for (pos = 0; pos < EXPAND_RATIO; pos++)
      {
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1)) + pos] +=
            ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos];
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1)) + pos] &=
            SecretInnerMask;

        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1) + 1) + pos] +=
            ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + (NbPlayers / 2) + i) + pos];
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1) + 1) + pos] &=
            SecretInnerMask;
      }
      for (; pos < 2 * EXPAND_RATIO; pos++)
      {
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1)) + pos] ^=
            ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos];
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1) + 1) + pos] ^=
            ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + (NbPlayers / 2) + i) + pos];
      }
    }
  }

  for (r = 0; r < NumRounds; r++)
  {
#pragma GCC unroll 4
#pragma GCC ivdep
    for (pos = 0; pos < EXPAND_RATIO; pos++)
    {
      Offsets[2 * EXPAND_RATIO * r + pos] =
          FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1)) + pos] +
          FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1) + 1) + pos];
      Offsets[2 * EXPAND_RATIO * r + pos] &= SecretInnerMask;
    }
    for (; pos < 2 * EXPAND_RATIO; pos++)
    {
      Offsets[2 * EXPAND_RATIO * r + pos] =
          FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1)) + pos] ^
          FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + LogPlayers - 1) + 1) + pos];
    }
  }
  for (int bitpos = 0; bitpos < LogPlayers - 1; bitpos++)
  {
    for (r = 0; r < NumRounds; r++)
    {
      for (pos = 0; pos < EXPAND_RATIO; pos++)
      {
        u64 opp;
        opp = -FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos];
        opp += CarryAddMask;
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos) + 1) + pos] = Offsets[2 * EXPAND_RATIO * r + pos] + opp;
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos) + 1) + pos] &= SecretInnerMask;
      }
      for (; pos < 2 * EXPAND_RATIO; pos++)
      {
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos) + 1) + pos] =
            FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos] ^
            Offsets[2 * EXPAND_RATIO * r + pos];
      }
    }
  }
  ClearUnused(NumRounds, FoldedExpTree, Offsets);
}

/*! \brief Hypercube folding of the expanded leaves
 */

// Need modular fold here where adequate

#ifdef CLASSICAL_FOLDING
void FastFoldTrees(int NumRounds, u64 *ExpandedTree, u64 *FoldedExpTree, u64 *Offsets)
#else
void FoldTrees(int NumRounds, u64 *ExpandedTree, u64 *FoldedExpTree, u64 *Offsets)
#endif
{
  int i, r, pos;

  // Folding [Left of each folded sharing]
  for (i = 0; i < 4 * EXPAND_RATIO * (NumRounds * LogPlayers); i++)
    FoldedExpTree[i] = 0;
  for (int bitpos = LogPlayers - 1; bitpos >= 0; bitpos--)
  {
    for (r = 0; r < NumRounds; r++)
    {
      for (i = 0; i < (1L << bitpos); i++)
      {
#pragma GCC unroll 4
#pragma GCC ivdep
        for (pos = 0; pos < EXPAND_RATIO; pos++)
        {
          FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos] +=
              ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos];
          ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos] +=
              ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i + (1L << bitpos)) + pos];
          FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos] &= SecretInnerMask;
          ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos] &= SecretInnerMask;
        }
        for (; pos < 2 * EXPAND_RATIO; pos++)
        {
          FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos] ^=
              ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos];
          ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i) + pos] ^=
              ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + i + (1L << bitpos)) + pos];
        }
      }
    }
  }
  // Finish the folding [Right of each folded sharing]
  // At this point ExpandedTree contains the full sum => Copy to Offsets
  for (r = 0; r < NumRounds; r++)
  {
    for (pos = 0; pos < 2 * EXPAND_RATIO; pos++)
    {
      Offsets[2 * EXPAND_RATIO * r + pos] = ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers) + pos];
    }
  }

  for (int bitpos = LogPlayers - 1; bitpos >= 0; bitpos--)
  {
    for (r = 0; r < NumRounds; r++)
    {
#pragma GCC unroll 4
#pragma GCC ivdep
      for (pos = 0; pos < EXPAND_RATIO; pos++)
      {
        u64 opp;
        opp = -FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos];
        opp += CarryAddMask;
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos) + 1) + pos] =
            (opp + Offsets[2 * EXPAND_RATIO * r + pos]) & SecretInnerMask;
      }
      for (; pos < 2 * EXPAND_RATIO; pos++)
      {
        FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos) + 1) + pos] =
            FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos)) + pos] ^
            Offsets[2 * EXPAND_RATIO * r + pos];
      }
    }
  }

  ClearUnused(NumRounds, FoldedExpTree, Offsets);
}

/*! \brief Recursive Fusion Sort calls
 */
void FusionRecursAsc(u64 *List, u64 listsize, u64 *Aux)
{
  u64 listsize1, i1, i2, i;

  if (listsize <= 1)
    return;

  listsize1 = (listsize + 1) / 2;
  FusionRecursAsc(Aux, listsize1, List);
  FusionRecursAsc(Aux + listsize1, listsize - listsize1, List + listsize1);
  for (i1 = 0, i2 = listsize1, i = 0; (i1 < listsize1) && (i2 < listsize); i++)
  {
    if (Aux[i1] <= Aux[i2])
    {
      List[i] = Aux[i1];
      i1++;
    }
    else
    {
      List[i] = Aux[i2];
      i2++;
    }
  }
  for (; i1 < listsize1; i1++, i++)
    List[i] = Aux[i1];
  for (; i2 < listsize; i2++, i++)
    List[i] = Aux[i2];
}

/*! \brief Recursive Fusion Sort entry point
 */
void ListFusionSortAsc(u64 *List, u64 listsize, u64 *Aux)
{
  u64 i;

  if (listsize <= 1)
    return;
  for (i = 0; i < listsize; i++)
    Aux[i] = List[i];
  FusionRecursAsc(List, listsize, Aux);
}

/*! \brief Create Prepermutation randomness by hashing
 */
void GetRawPermCoefficients(u64 *PermTree, aessubkeytype *key_schedules, int lvl_offset)
{
  u64 CurrentSize;
  for (CurrentSize = 2; CurrentSize < PERM_LEAVES; CurrentSize <<= 1)
  {
    aes128_dbl_enc(PermTree + 4 * (PERM_LEAVES - CurrentSize),
                   PermTree + 4 * (PERM_LEAVES - 2 * CurrentSize),
                   CurrentSize, key_schedules + 22 * lvl_offset,
                   key_schedules + 22 * lvl_offset + 11);
  }
}

/*! \brief Make Permutations
 */
void MakePermutations(u64 *Permuts, u64 *PermTree, aessubkeytype *key_schedules, int lvl_offset)
{
  u64 r, i;
  u64 permmask = 0xffffffffffffff00UL;

  GetRawPermCoefficients(PermTree, key_schedules, lvl_offset);

  for (r = 0; r < NbRounds; r++)
  {
    for (i = 0; i < NumWindows; i++)
    {
      Permuts[r * NumWindows + i] = (PermTree[r * NumWindows + i] & permmask) | i;
    }
    ListFusionSortAsc(Permuts + r * NumWindows, NumWindows, Permuts + NbRounds * NumWindows);
    for (i = 0; i < NumWindows; i++)
    {
      Permuts[r * NumWindows + i] &= 0xff;
    }
  }
}

/*! \brief Copy a byte to the signature output buffer
 */
void CharsToSignature(uchar *buff, int size, uchar *CompSig, int *sigsize)
{
  int i;
  for (i = 0; i < size; i++)
    CompSig[(*sigsize) + i] = buff[i];
  *sigsize += size;
}

/*! \brief Copy a 64-bit integer to the signature output buffer
 */
void U64ToSignature(u64 *buff, int size, uchar *CompSig, int *sigsize)
{
  int cu64, cbyte, offset;

  offset = *sigsize;
  for (cu64 = 0; cu64 < size; cu64++)
  {
    for (cbyte = 0; cbyte < 8; cbyte++)
    {
      CompSig[offset] = (buff[cu64] >> (8 * cbyte)) & 0xff;
      offset++;
    }
  }
  *sigsize = offset;
}

/*! \brief Copy a 64-bit integer (possibly partial) to the signature output buffer
 */
void U64ChunksToSignature(u64 *buff, int size, uchar *CompSig, int *sigsize)
{
  int cu64, cbyte, offset, count;

  offset = *sigsize;
  cu64 = 0;
  cbyte = 0;
  for (count = 0; count < size; count++)
  {
    CompSig[offset] = (buff[cu64] >> (8 * cbyte)) & 0xff;
    offset++;
    cbyte++;
    if (cbyte == 8)
    {
      cu64++;
      cbyte = 0;
    }
  }
  *sigsize = offset;
}

/*! \brief Read a byte from the signature output buffer
 */
void CharsFromSignature(uchar *buff, int size, uchar *CompSig, int *sigsize)
{
  int i;
  for (i = 0; i < size; i++)
    buff[i] = CompSig[(*sigsize) + i];
  *sigsize += size;
}

/*! \brief Read a 64-bit integer from the signature output buffer
 */
void U64FromSignature(u64 *buff, int size, uchar *CompSig, int *sigsize)
{
  int cu64, cbyte, offset;

  offset = *sigsize;
  for (cu64 = 0; cu64 < size; cu64++)
  {
    buff[cu64] = 0;
    for (cbyte = 0; cbyte < 8; cbyte++)
    {
      buff[cu64] |= ((u64)CompSig[offset]) << (8 * cbyte);
      offset++;
    }
  }
  *sigsize = offset;
}

/*! \brief Read a 64-bit integer (possible partial) from the signature output buffer
 */
void U64ChunksFromSignature(u64 *buff, int size, uchar *CompSig, int *sigsize)
{
  int cu64, cbyte, offset, count;

  offset = *sigsize;
  cu64 = 0;
  cbyte = 0;
  buff[cu64] = 0;
  for (count = 0; count < size; count++)
  {
    buff[cu64] |= ((u64)CompSig[offset]) << (8 * cbyte);
    offset++;
    cbyte++;
    if (cbyte == 8)
    {
      cu64++;
      buff[cu64] = 0;
      cbyte = 0;
    }
  }
  *sigsize = offset;
}

#ifdef RSD_EVALUATE
void OpenPermuted(u64 *FoldedExpTree, u64 *Offsets, u64 *Permuts, u64 *PermutedOpenings)
{
  int r, i, j;
  u64 ModsumShares[2 * SecretRepinU128];

  for (r = 0; r < NbRounds; r++)
  {
    for (j = 0; j < 2 * SecretRepinU128; j++)
    {
      PermutedOpenings[r * 2 * SecretRepinU128 + j] = Offsets[2 * EXPAND_RATIO * r + j];
    }
    for (i = 0; i < 2; i++)
    {
      for (j = 0; j < 2 * SecretRepinU128; j++)
      {
        ModsumShares[j] = 0;
      }
      for (j = 0; j < NumWindows; j++)
      {
        int inpos;
        u64 val;
        inpos = Permuts[r * NumWindows + j];
        val = (FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + 2 * SecretRepinU128 + (inpos / 16)] >> ((inpos & 15) * 4)) & 0x7;
        ModsumShares[j / 16] |= val << ((j & 15) * 4);
      }
      for (j = 0; j < 2 * SecretRepinU128; j++)
      {
        PermutedOpenings[r * 2 * SecretRepinU128 + j] += ModsumShares[j];
        PermutedOpenings[r * 2 * SecretRepinU128 + j] &= SecretInnerMask;
        PermutedOpenings[r * 2 * SecretRepinU128 + j] += FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + j];
        PermutedOpenings[r * 2 * SecretRepinU128 + j] &= SecretInnerMask;
      }
    }
  }
}

void CompleteParities(u64 *FoldedExpTree, u64 *Offsets)
{
  int r, i, j;
  u64 byteparity, tmp;
  for (r = 0; r < NbRounds; r++)
  {
    for (j = 0; j < EXPAND_RATIO; j++)
    {
      tmp = Offsets[2 * EXPAND_RATIO * r + EXPAND_RATIO + j];
      tmp = (tmp ^ (tmp >> 4)) & CompressHalf;
      tmp = (tmp ^ (tmp >> 2)) & CompressQuat;
      tmp = (tmp ^ (tmp >> 1)) & CompressLow;
      Offsets[2 * EXPAND_RATIO * r + EXPAND_RATIO + j] ^= CompressHigh ^ (tmp << 7);
    }
    for (i = 0; i < 2 * LogPlayers; i++)
    {
      for (j = 0; j < EXPAND_RATIO; j++)
      {
        tmp = FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + EXPAND_RATIO + j];
        tmp = (tmp ^ (tmp >> 4)) & CompressHalf;
        tmp = (tmp ^ (tmp >> 2)) & CompressQuat;
        tmp = (tmp ^ (tmp >> 1)) & CompressLow;
        FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + EXPAND_RATIO + j] ^= (tmp << 7);
      }
    }
  }
}

/*! \brief Evaluate and hash players' views
 */
void DoRSDevaluate(SHA256_CTX *hash_ptr, u64 *Synd, u64 *PubMat, u64 *FoldedExpTree,
                   u64 *Offsets, u64 *Permuts, u64 *PermutedOpenings)
{
  int r, i, j;

#ifdef COMPRESSION
  CompleteParities(FoldedExpTree, Offsets);
#endif
  OpenPermuted(FoldedExpTree, Offsets, Permuts, PermutedOpenings);

  for (r = 0; r < NbRounds; r++)
  {
    u64 SyndShares[2 * CodeRepinU128];
    u64 SyndSharesSum[2 * CodeRepinU128];
    for (i = 0; i < 2 * LogPlayers; i++)
    {
      u64 XShares[EXPAND_RATIO];
      u64 ModsumShares[2 * SecretRepinU128];

      // Commit x-Pi(z) shares [Use addition// see z/u generate]
      for (j = 0; j < 2 * SecretRepinU128; j++)
      {
        ModsumShares[j] = 0;
      }
      for (j = 0; j < NumWindows; j++)
      {
        int inpos;
        u64 val;
        inpos = Permuts[r * NumWindows + j];
        val = (FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + 2 * SecretRepinU128 + (inpos / 16)] >> ((inpos & 15) * 4)) & 0x7;
        ModsumShares[j / 16] |= val << ((j & 15) * 4);
      }
      for (j = 0; j < 2 * SecretRepinU128; j++)
      {
        ModsumShares[j] += FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + j];
        ModsumShares[j] &= SecretInnerMask;
      }
      SHA256Update(hash_ptr, (uchar *)ModsumShares, 2 * SecretRepinU128 * sizeof(u64));
      //    printf("Modshares[%d,%d]:%x\n",r,i,hash_ptr->state.st32[0]);

#ifdef OPTIM_SCALARS
      if ((i > 1) && (i & 1))
      {
        SHA256Update(hash_ptr, (uchar *)SyndShares, 2 * CodeRepinU128 * sizeof(u64));
        //	printf("Syndshares[%d,%d]:%x\n",r,i,hash_ptr->state.st32[0]);
        continue; // Compute fewer scalar products
      }
#endif

      for (j = 0; j < EXPAND_RATIO; j++)
      {
        u64 bitpattern = 0;
        for (int w = 0; w < 8; w++)
        {
          u64 tmp;
          int inpos, permpos, shift;
          inpos = j * 8 + w;
          if (inpos < NumWindows)
          {
            permpos = Permuts[r * NumWindows + inpos];
            tmp = FoldedExpTree[2 * EXPAND_RATIO * (2 * r * LogPlayers + i) + EXPAND_RATIO + (permpos / 8)];
            tmp = (tmp >> (8 * (permpos & 7))) & 0xff;
            shift = (PermutedOpenings[r * 2 * SecretRepinU128 + (inpos / 16)] >> ((inpos & 15) * 4)) & 7;
            tmp = ((tmp << shift) | (tmp >> (8 - shift))) & 0xff;
            bitpattern |= tmp << (8 * w);
          }
        }
        XShares[j] = bitpattern;
      }

      // Multiply Transformed share
      DenseMatMul(SyndShares, XShares, PubMat);
      // Commit Syndrome shares
      SHA256Update(hash_ptr, (uchar *)SyndShares, 2 * CodeRepinU128 * sizeof(u64));
      //  printf("Syndshares[%d,%d]:%x\n",r,i,hash_ptr->state.st32[0]);

#ifdef OPTIM_SCALARS
      if (i == 0)
      {
        for (j = 0; j < 2 * CodeRepinU128; j++)
        {
          SyndSharesSum[j] = SyndShares[j];
        }
      }
      else if (i == 1)
      {
        for (j = 0; j < 2 * CodeRepinU128; j++)
        {
          SyndSharesSum[j] ^= SyndShares[j];
        }
      }
      else
      {
        for (j = 0; j < 2 * CodeRepinU128; j++)
        {
          SyndShares[j] ^= SyndSharesSum[j];
        }
      }
#endif
    }
  }
}

/*! \brief Check the hashed players' views
 */
void DoRSDevaluateCheck(SHA256_CTX *hash_ptr, int *Queries, u64 *Synd, u64 *PubMat, u64 *FoldedExpTree,
                        u64 *Offsets, u64 *Permuts, u64 *PermutedOpenings)
{
  int r, i, j;

#ifdef COMPRESSION
  CompleteParities(FoldedExpTree, Offsets);
#endif

  for (r = 0; r < NbRounds; r++)
  {
    u64 Modsum[2 * SecretRepinU128];
    u64 SyndDiff[2 * CodeRepinU128];
    u64 XShares[EXPAND_RATIO];

    for (j = 0; j < EXPAND_RATIO; j++)
    {
      u64 bitpattern = 0;
      for (int w = 0; w < 8; w++)
      {
        u64 tmp;
        int inpos, permpos, shift;
        inpos = j * 8 + w;
        if (inpos < NumWindows)
        {
          permpos = Permuts[r * NumWindows + inpos];
          tmp = Offsets[2 * EXPAND_RATIO * r + EXPAND_RATIO + (permpos / 8)];
          tmp = (tmp >> (8 * (permpos & 7))) & 0xff;
          shift = (PermutedOpenings[r * 2 * SecretRepinU128 + (inpos / 16)] >> ((inpos & 15) * 4)) & 7;
          tmp = ((tmp << shift) | (tmp >> (8 - shift))) & 0xff;
          bitpattern |= tmp << (8 * w);
        }
      }
      XShares[j] = bitpattern;
    }

    for (j = 0; j < 2 * CodeRepinU128; j++)
      SyndDiff[j] = Synd[j];
    DenseMatMulAdd(SyndDiff, XShares, PubMat);

    for (j = 0; j < 2 * SecretRepinU128; j++)
    {
      Modsum[j] = (-Offsets[2 * EXPAND_RATIO * r + j]) + CarryAddMask;
      Modsum[j] &= SecretInnerMask;
      Modsum[j] += PermutedOpenings[r * 2 * SecretRepinU128 + j];
      Modsum[j] &= SecretInnerMask;
    }

    for (int bitpos = 0; bitpos < LogPlayers; bitpos++)
    {
      int good, bad, pos;
      u64 SyndShares[4 * CodeRepinU128];
      u64 ModsumShares[4 * SecretRepinU128];

      bad = ((Queries[r] >> bitpos) & 1);
      good = bad ^ 1;

      for (j = 0; j < 4 * SecretRepinU128; j++)
      {
        ModsumShares[j] = 0;
      }
      for (j = 0; j < NumWindows; j++)
      {
        int inpos;
        u64 val;
        inpos = Permuts[r * NumWindows + j];
        val = (FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos) + good) + 2 * SecretRepinU128 + (inpos / 16)] >> ((inpos & 15) * 4)) & 0x7;
        ModsumShares[2 * SecretRepinU128 * good + (j / 16)] |= val << ((j & 15) * 4);
      }

      // Commit x-Pi(z) shares [Use addition// see z/u generate]
      for (j = 0; j < 2 * SecretRepinU128; j++)
      {
        ModsumShares[2 * SecretRepinU128 * good + j] += FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos) + good) + j];
        ModsumShares[2 * SecretRepinU128 * good + j] &= SecretInnerMask;
        ModsumShares[2 * SecretRepinU128 * bad + j] = (-ModsumShares[2 * SecretRepinU128 * good + j] + CarryAddMask) & SecretInnerMask;
        ModsumShares[2 * SecretRepinU128 * bad + j] += Modsum[j];
        ModsumShares[2 * SecretRepinU128 * bad + j] &= SecretInnerMask;
      }

      for (j = 0; j < EXPAND_RATIO; j++)
      {
        u64 bitpattern = 0;
        for (int w = 0; w < 8; w++)
        {
          u64 tmp;
          int inpos, permpos, shift;
          inpos = j * 8 + w;
          if (inpos < NumWindows)
          {
            permpos = Permuts[r * NumWindows + inpos];
            tmp = FoldedExpTree[2 * EXPAND_RATIO * (2 * (r * LogPlayers + bitpos) + good) + EXPAND_RATIO + (permpos / 8)];
            tmp = (tmp >> (8 * (permpos & 7))) & 0xff;
            shift = (PermutedOpenings[r * 2 * SecretRepinU128 + (inpos / 16)] >> ((inpos & 15) * 4)) & 7;
            tmp = ((tmp << shift) | (tmp >> (8 - shift))) & 0xff;
            bitpattern |= tmp << (8 * w);
          }
        }
        XShares[j] = bitpattern;
      }

      // Multiply Transformed share and get bad one
      DenseMatMul(SyndShares + 2 * CodeRepinU128 * good, XShares, PubMat);
      for (j = 0; j < 2 * CodeRepinU128; j++)
      {
        SyndShares[2 * CodeRepinU128 * bad + j] = SyndShares[2 * CodeRepinU128 * good + j] ^ SyndDiff[j];
      }

      // Redo commitments in same order
      SHA256Update(hash_ptr, (uchar *)ModsumShares, 2 * SecretRepinU128 * sizeof(u64));
      //  printf("Modshares[%d,%d]:%x\n",r,2*bitpos,hash_ptr->state.st32[0]);
      SHA256Update(hash_ptr, (uchar *)SyndShares, 2 * CodeRepinU128 * sizeof(u64));
      //    printf("Syndshares[%d,%d]:%x\n",r,2*bitpos,hash_ptr->state.st32[0]);
      SHA256Update(hash_ptr, (uchar *)(ModsumShares + 2 * SecretRepinU128), 2 * SecretRepinU128 * sizeof(u64));
      //   printf("Modshares[%d,%d]:%x\n",r,2*bitpos+1,hash_ptr->state.st32[0]);
      SHA256Update(hash_ptr, (uchar *)(SyndShares + 2 * CodeRepinU128), 2 * CodeRepinU128 * sizeof(u64));
      //   printf("Syndshares[%d,%d]:%x\n",r,2*bitpos+1,hash_ptr->state.st32[0]);
    }
  }
}
#endif

/*! \brief Main routine to compute a RSD signature
 */
void DoRSDSignature(uchar *messhash, uchar *randomness, u64 *SecX,
                    u64 *Synd, u64 *PubMat, uchar *CompSig, int *sigsize,
                    u64 *workMemory)
// Inputs: hash of message to sign, random seed, Secret key
//         Expanded PK, SKXPK scalar products, Signature buffer,
//         Signature size (returned)
//         Preallocsated working memory
{
  SHA256_CTX hash_ctx, hash_ctx_cpy;
  uchar inner_random[32 * NB_INNER_RANDOM];
  u64 AES_keys[2 * NB_AES_KEYS];
  aessubkeytype key_schedules[11 * NB_AES_KEYS];
  u64 *FullTree, *ExpandedTree, *FoldedExpTree, *PermTree;
  u64 *Offsets, *Permuts, *PermutedOpenings;
  //  u64 *OffsetsCmp, *FoldedExpTreeCmp;
  u64 *TmpAllocPtr, CurrentSize;

  int r, i, lvl;
  int Queries[NbRounds];

  // Create pointers
  TmpAllocPtr = workMemory;
  Offsets = AssignNewU64Ptr(TmpAllocPtr, 2 * EXPAND_RATIO * NbRounds);
  //  OffsetsCmp=AssignNewU64Ptr(TmpAllocPtr,2*EXPAND_RATIO*NbRounds);
  FullTree = AssignNewU64Ptr(TmpAllocPtr, 4 * LEAVES_SIZE);
  ExpandedTree = AssignNewU64Ptr(TmpAllocPtr, 2 * EXPAND_RATIO * LEAVES_SIZE);
  FoldedExpTree = AssignNewU64Ptr(TmpAllocPtr, 4 * EXPAND_RATIO * NbRounds * LogPlayers);
  // FoldedExpTreeCmp=AssignNewU64Ptr(TmpAllocPtr,4*EXPAND_RATIO*NbRounds*LogPlayers);
  PermTree = AssignNewU64Ptr(TmpAllocPtr, 4 * PERM_LEAVES);
  Permuts = AssignNewU64Ptr(TmpAllocPtr, (NbRounds + 1) * NumWindows);
  PermutedOpenings = AssignNewU64Ptr(TmpAllocPtr, 2 * SecretRepinU128 * NbRounds);

  InitKeyMaterial(messhash, randomness,
                  AES_keys, key_schedules, NB_AES_KEYS,
                  inner_random, NB_INNER_RANDOM);

  InitTrees(NbRounds, SecX, 0, inner_random, FullTree);

  GenerateTrees(NbRounds, 0, key_schedules,
                FullTree, ExpandedTree, NULL);

  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, messhash, 32);
  SHA256Update(&hash_ctx, (uchar *)(AES_keys), 2 * sizeof(u64)); // K0

  //  SlowFoldTrees(NbRounds, ExpandedTree, FoldedExpTreeCmp, OffsetsCmp);
  FoldTrees(NbRounds, ExpandedTree, FoldedExpTree, Offsets);

  // Modular offset for sharing of X
  for (r = 0; r < NbRounds; r++)
  {
    for (i = 0; i < 2 * SecretRepinU128; i++)
    {
      u64 opp;
      opp = -Offsets[2 * EXPAND_RATIO * r + i];
      opp += CarryAddMask;
      Offsets[2 * EXPAND_RATIO * r + i] = (SecX[i] + opp) & SecretInnerMask;
    }
  }

  // Make random full sharings conform to random modular sharings
  for (r = 0; r < NbRounds; r++)
  {
    for (i = 0; i < EXPAND_RATIO; i++)
    {
      int pos;
      u64 bitpattern = 0;
      for (pos = 0; pos < 8; pos++)
      {
        int modshare;
        modshare = (Offsets[2 * EXPAND_RATIO * r + 2 * SecretRepinU128 + (i / 2)] >> (((i & 1) * 8 + pos) * 4)) & 7;
        modshare = (8 - modshare) & 7;
        if (8 * i + pos < NumWindows)
          bitpattern |= 1UL << (8 * pos + modshare);
      }
      Offsets[2 * EXPAND_RATIO * r + EXPAND_RATIO + i] ^= bitpattern;
#ifdef COMPRESSION
      Offsets[2 * EXPAND_RATIO * r + EXPAND_RATIO + i] &= CompressMask;
#endif
    }
  }

  for (r = 0; r < NbRounds; r++)
  {
    SHA256Update(&hash_ctx, (uchar *)(Offsets + 2 * EXPAND_RATIO * r), 2 * SecretRepinU128 * sizeof(u64));
  }

  for (r = 0; r < NbRounds; r++)
  {
    SHA256Update(&hash_ctx, (uchar *)(Offsets + 2 * EXPAND_RATIO * r + EXPAND_RATIO), EXPAND_RATIO * sizeof(u64));
  }
  hash_ctx_cpy = hash_ctx;

  SHA256Final((uchar *)(PermTree + 4 * (PERM_LEAVES - 2)), &hash_ctx);
  MakePermutations(Permuts, PermTree, key_schedules, LogPlayers + 1);

  // Change parameters and content [commit syndrome shares]
  DoRSDevaluate(&hash_ctx_cpy, Synd, PubMat, FoldedExpTree, Offsets, Permuts, PermutedOpenings);

  // Hash to Signature and Get Queries
  SHA256Final(CompSig, &hash_ctx_cpy);
  *sigsize = 32;
  for (r = 0; r < NbRounds; r++)
  {
#ifdef L8
    Queries[r] = CompSig[r];
#else
    Queries[r] = ((((int)CompSig[2 * r]) << 8) + CompSig[2 * r + 1]) & ((1 << LogPlayers) - 1);
#endif
  }
  // Add K0 to signature
  U64ToSignature(AES_keys, 2, CompSig, sigsize);

  for (r = 0; r < NbRounds; r++)
  {
    // Punctured key(s) to signature
    CurrentSize = NbRounds; // Level with NbRounds nodes is last unread
    for (lvl = 0; lvl < LogPlayers; lvl++)
    {
      u64 *ReadTree;
      int ReadPos, ZeroPos;
      CurrentSize <<= 1;
      ReadTree = FullTree + 4 * (LEAVES_SIZE - CurrentSize);
      ReadPos = ((Queries[r] >> (LogPlayers - 1 - lvl)) ^ 1) + r * (CurrentSize / NbRounds);
      ZeroPos = ((Queries[r] >> (LogPlayers - 1 - lvl))) + r * (CurrentSize / NbRounds);
      U64ToSignature(ReadTree + 2 * ReadPos, 2, CompSig, sigsize);
#ifdef SANITIZE
      ReadTree[2 * ZeroPos] = 0;
      ReadTree[2 * ZeroPos + 1] = 0;
#endif
    }
#ifdef SANITIZE
    for (lvl = 0; lvl < LogPlayers; lvl++)
    {
      int foldedquery;
      foldedquery = Queries[r] & ((2 << lvl) - 1);
      for (i = 0; i < 2 * EXPAND_RATIO; i++)
        ExpandedTree[2 * EXPAND_RATIO * (r * NbPlayers + foldedquery) + i] = 0;
    }
    for (i = 0; i < 4 * EXPAND_RATIO * NbRounds * LogPlayers; i++)
      FoldedExpTree[i] = 0;
#endif
  }

  // Global Elements to signature
#ifndef COMPRESSION
  for (r = 0; r < NbRounds; r++)
  {
    U64ChunksToSignature(Offsets + 2 * EXPAND_RATIO * r + EXPAND_RATIO, NumWindows, CompSig, sigsize);
  }
  {
    uchar bitbuffer[(NbRounds * NumWindows * 6 + 7) / 8];
    int count;
    for (count = 0; count < (NbRounds * NumWindows * 6 + 7) / 8; count++)
      bitbuffer[count] = 0;
    count = 0;
    for (r = 0; r < NbRounds; r++)
    {
      for (int w = 0; w < NumWindows; w++)
      {
        int wbyte;
        wbyte = (Offsets[2 * EXPAND_RATIO * r + (w / 16)] >> (4 * (w & 15))) & 0x7;
        wbyte = (wbyte << 3) | ((PermutedOpenings[2 * SecretRepinU128 * r + (w / 16)] >> (4 * (w & 15))) & 0x7);
        for (int bb = 0; bb < 6; bb++)
        {
          bitbuffer[count / 8] ^= ((wbyte >> bb) & 1) << (count & 7);
          count++;
        }
      }
    }
    CharsToSignature(bitbuffer, (count + 7) / 8, CompSig, sigsize);
  }
#else
  {
    uchar bitbuffer[(NbRounds * NumWindows * (6 + 7) + 7) / 8];
    int count;
    for (count = 0; count < (NbRounds * NumWindows * (6 + 7) + 7) / 8; count++)
      bitbuffer[count] = 0;
    count = 0;
    for (r = 0; r < NbRounds; r++)
    {
      for (int w = 0; w < NumWindows; w++)
      {
        int wbyte;
        wbyte = (Offsets[2 * EXPAND_RATIO * r + EXPAND_RATIO + (w / 8)] >> (8 * (w & 7))) & 0x7f;
        for (int bb = 0; bb < 7; bb++)
        {
          bitbuffer[count / 8] ^= ((wbyte >> bb) & 1) << (count & 7);
          count++;
        }
      }
    }
    for (r = 0; r < NbRounds; r++)
    {
      for (int w = 0; w < NumWindows; w++)
      {
        int wbyte;
        wbyte = (Offsets[2 * EXPAND_RATIO * r + (w / 16)] >> (4 * (w & 15))) & 0x7;
        wbyte = (wbyte << 3) | ((PermutedOpenings[2 * SecretRepinU128 * r + (w / 16)] >> (4 * (w & 15))) & 0x7);
        for (int bb = 0; bb < 6; bb++)
        {
          bitbuffer[count / 8] ^= ((wbyte >> bb) & 1) << (count & 7);
          count++;
        }
      }
    }
    CharsToSignature(bitbuffer, (count + 7) / 8, CompSig, sigsize);
  }
#endif
}

/*! \brief Main routine to verify a RSD signature
 */
int CheckRSDSignature(uchar *messhash, u64 *Synd, u64 *PubMat, uchar *CompSig,
                      u64 *workMemory)

{
  SHA256_CTX hash_ctx, hash_ctx_cpy;
  u64 AES_keys[2 * NB_AES_KEYS];
  aessubkeytype key_schedules[11 * NB_AES_KEYS];
  u64 *FullTree, *ExpandedTree, *FoldedExpTree, *PermTree;
  u64 *Offsets, *Permuts, *PermutedOpenings;
  u64 *TmpAllocPtr, CurrentSize;
  uchar Tmp[32];

  int r, i, lvl, ok;
  int sigoffset;
  int Queries[NbRounds];

  TmpAllocPtr = workMemory;
  Offsets = AssignNewU64Ptr(TmpAllocPtr, 2 * EXPAND_RATIO * NbRounds);
  FullTree = AssignNewU64Ptr(TmpAllocPtr, 4 * LEAVES_SIZE);
  ExpandedTree = AssignNewU64Ptr(TmpAllocPtr, 2 * EXPAND_RATIO * LEAVES_SIZE);
  FoldedExpTree = AssignNewU64Ptr(TmpAllocPtr, 4 * EXPAND_RATIO * NbRounds * LogPlayers);
  PermTree = AssignNewU64Ptr(TmpAllocPtr, 4 * PERM_LEAVES);
  Permuts = AssignNewU64Ptr(TmpAllocPtr, (NbRounds + 1) * NumWindows);
  PermutedOpenings = AssignNewU64Ptr(TmpAllocPtr, 2 * SecretRepinU128 * NbRounds);

  for (r = 0; r < NbRounds; r++)
  {
#ifdef L8
    Queries[r] = CompSig[r];
#else
    Queries[r] = ((((int)CompSig[2 * r]) << 8) + CompSig[2 * r + 1]) & ((1 << LogPlayers) - 1);
#endif
  }
  sigoffset = 32;
  // Get K0 From signature
  U64FromSignature(AES_keys, 2, CompSig, &sigoffset);
  FillAESKey(AES_keys, key_schedules, NB_AES_KEYS);

  for (r = 0; r < NbRounds; r++)
  {
    // Punctured key(s) to signature
    CurrentSize = NbRounds; // Level with NbRounds nodes is last unread
    for (lvl = 0; lvl < LogPlayers; lvl++)
    {
      u64 *ReadTree;
      int ReadPos;
      CurrentSize <<= 1;
      ReadTree = FullTree + 4 * (LEAVES_SIZE - CurrentSize);
      ReadPos = ((Queries[r] >> (LogPlayers - 1 - lvl)) ^ 1) + r * (CurrentSize / NbRounds);
      U64FromSignature(ReadTree + 2 * ReadPos, 2, CompSig, &sigoffset);
    }
  }

  // Reconstruct tree from path already in it
  GenerateTrees(NbRounds, 0,
                key_schedules,
                FullTree, ExpandedTree, Queries);

  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, messhash, 32);
  SHA256Update(&hash_ctx, (uchar *)(AES_keys), 2 * sizeof(u64)); // K0

  FoldTrees(NbRounds, ExpandedTree, FoldedExpTree, Offsets);

  for (i = 0; i < 2 * EXPAND_RATIO * NbRounds; i++)
    Offsets[i] = 0;
  for (i = 0; i < 2 * SecretRepinU128 * NbRounds; i++)
    PermutedOpenings[i] = 0;

    // Global Elements from signature
#ifndef COMPRESSION
  for (r = 0; r < NbRounds; r++)
  {
    U64ChunksFromSignature(Offsets + 2 * EXPAND_RATIO * r + EXPAND_RATIO, NumWindows, CompSig, &sigoffset);
  }
  {
    uchar bitbuffer[(NbRounds * NumWindows * 6 + 7) / 8];
    int count;
    CharsFromSignature(bitbuffer, (NbRounds * NumWindows * 6 + 7) / 8, CompSig, &sigoffset);
    count = 0;
    for (r = 0; r < NbRounds; r++)
    {
      for (int w = 0; w < NumWindows; w++)
      {
        u64 wbyte;
        wbyte = 0;
        for (int bb = 0; bb < 6; bb++)
        {
          wbyte |= ((bitbuffer[count / 8] >> (count & 7)) & 1) << bb;
          count++;
        }
        Offsets[2 * EXPAND_RATIO * r + (w / 16)] ^= (wbyte >> 3) << (4 * (w & 15));
        PermutedOpenings[2 * SecretRepinU128 * r + (w / 16)] ^= (wbyte & 0x7UL) << (4 * (w & 15));
      }
    }
  }
#else
  {
    uchar bitbuffer[(NbRounds * NumWindows * (6 + 7) + 7) / 8];
    int count;
    CharsFromSignature(bitbuffer, (NbRounds * NumWindows * 7 + 7) / 8, CompSig, &sigoffset);
    count = 0;
    for (r = 0; r < NbRounds; r++)
    {
      for (int w = 0; w < NumWindows; w++)
      {
        u64 wbyte;
        wbyte = 0;
        for (int bb = 0; bb < 7; bb++)
        {
          wbyte |= ((bitbuffer[count / 8] >> (count & 7)) & 1) << bb;
          count++;
        }
        Offsets[2 * EXPAND_RATIO * r + EXPAND_RATIO + (w / 8)] ^= (wbyte) << (8 * (w & 7));
      }
    }
    for (r = 0; r < NbRounds; r++)
    {
      for (int w = 0; w < NumWindows; w++)
      {
        u64 wbyte;
        wbyte = 0;
        for (int bb = 0; bb < 6; bb++)
        {
          wbyte |= ((bitbuffer[count / 8] >> (count & 7)) & 1) << bb;
          count++;
        }
        Offsets[2 * EXPAND_RATIO * r + (w / 16)] ^= (wbyte >> 3) << (4 * (w & 15));
        PermutedOpenings[2 * SecretRepinU128 * r + (w / 16)] ^= (wbyte & 0x7UL) << (4 * (w & 15));
      }
    }
  }
#endif

  for (r = 0; r < NbRounds; r++)
  {
    SHA256Update(&hash_ctx, (uchar *)(Offsets + 2 * EXPAND_RATIO * r), 2 * SecretRepinU128 * sizeof(u64));
  }

  for (r = 0; r < NbRounds; r++)
  {
#ifdef COMPRESSION
    for (int pos = EXPAND_RATIO; pos < 2 * EXPAND_RATIO; pos++)
      Offsets[2 * EXPAND_RATIO * r + pos] &= CompressMask;
#endif
    SHA256Update(&hash_ctx, (uchar *)(Offsets + 2 * EXPAND_RATIO * r + EXPAND_RATIO), EXPAND_RATIO * sizeof(u64));
  }
  hash_ctx_cpy = hash_ctx;
  SHA256Final((uchar *)(PermTree + 4 * (PERM_LEAVES - 2)), &hash_ctx);
  MakePermutations(Permuts, PermTree, key_schedules, LogPlayers + 1);

  DoRSDevaluateCheck(&hash_ctx_cpy, Queries, Synd, PubMat, FoldedExpTree, Offsets, Permuts, PermutedOpenings);

  ok = 1;
  SHA256Final(Tmp, &hash_ctx_cpy);
  for (i = 0; i < 32; i++)
  {
    if (Tmp[i] != CompSig[i])
      ok = 0;
  }
  return ok;
}

/*! \brief Main program. Test and time the signature scheme (nbloop repetitions)
 */
int main(int argc, char **argv)
{

  SHA256_CTX hash_ctx;
  u64 Seed[4];
  u64 SecX[2 * SecretRepinU128];
  u64 Synd[2 * CodeRepinU128];
  u64 PubMat[2 * CodeLength * CodeRepinU128];
  u64 PubSeed[2];
  int workMemSize;
  u64 *workMem;

  uchar messhash[32];
  uchar randomness[32];
  uchar CompSig[10000];
  int loop;
  int accum;
  int sigsize;

  int nbloop = 100;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define B_KEY_GENERATION 0
#define B_SIGN_ALGO 1
#define B_VERIFY_ALGO 2
#define NUMBER_OF_ALGO_BENCHES 3
#endif

  btimer_t timers_algos[NUMBER_OF_ALGO_BENCHES];
  double std_timer[NUMBER_OF_ALGO_BENCHES];

  double timer_pow2[NUMBER_OF_ALGO_BENCHES];
  for (int j = 0; j < NUMBER_OF_ALGO_BENCHES; j++)
  {
    btimer_init(&timers_algos[j]);
    timer_pow2[j] = 0;
  }

  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, (uchar *)"secret seed", 11);
  SHA256Final((uchar *)Seed, &hash_ctx);

  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, (uchar *)"Message", 7);
  SHA256Final(messhash, &hash_ctx);

  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, (uchar *)"BlaBlaRandom", 12);
  SHA256Final(randomness, &hash_ctx);

  workMemSize = 0;
  workMemSize += 2 * EXPAND_RATIO * NbRounds;
  // workMemSize+=2*EXPAND_RATIO*NbRounds;  // Twice to compare the two foldings
  workMemSize += 4 * LEAVES_SIZE;
  workMemSize += 2 * EXPAND_RATIO * LEAVES_SIZE;
  workMemSize += 4 * EXPAND_RATIO * NbRounds * LogPlayers;
  // workMemSize+=4*EXPAND_RATIO*NbRounds*LogPlayers; // Twice to compare the two foldings
  workMemSize += 4 * PERM_LEAVES;
  workMemSize += (NbRounds + 1) * NumWindows;
  workMemSize += 2 * SecretRepinU128 * NbRounds;

  workMem = (u64 *)malloc(workMemSize * sizeof(u64));

  printf("Work memory size is:%ld\n", workMemSize * sizeof(u64));

  for (loop = 0; loop < nbloop; loop++)
  {
    btimer_start(&timers_algos[B_KEY_GENERATION]);
    CreateKeyFromSeed(Seed, SecX, Synd, PubMat, PubSeed);
    btimer_end(&timers_algos[B_KEY_GENERATION]);
    btimer_count(&timers_algos[B_KEY_GENERATION]);
    timer_pow2[B_KEY_GENERATION] += pow(btimer_diff(&timers_algos[B_KEY_GENERATION]), 2) / nbloop;
  }

  accum = 1;
  for (loop = 0; loop < nbloop; loop++)
  {
    int ok;
    btimer_start(&timers_algos[B_SIGN_ALGO]);
    ok = 0;
    DoRSDSignature(messhash, randomness, SecX, Synd, PubMat,
                   CompSig, &sigsize, workMem);
    btimer_end(&timers_algos[B_SIGN_ALGO]);
    btimer_count(&timers_algos[B_SIGN_ALGO]);
    timer_pow2[B_SIGN_ALGO] += pow(btimer_diff(&timers_algos[B_SIGN_ALGO]), 2) / nbloop;
    randomness[0]++;
    btimer_start(&timers_algos[B_VERIFY_ALGO]);
    accum &= CheckRSDSignature(messhash, Synd, PubMat, CompSig, workMem);
    btimer_end(&timers_algos[B_VERIFY_ALGO]);
    btimer_count(&timers_algos[B_VERIFY_ALGO]);
    timer_pow2[B_VERIFY_ALGO] += pow(btimer_diff(&timers_algos[B_VERIFY_ALGO]), 2) / nbloop;
  }
  printf("Compressed size is:%d\n", sigsize);
  printf("Return cor: %d\n", accum);

  std_timer[B_KEY_GENERATION] = sqrt(timer_pow2[B_KEY_GENERATION] - pow(btimer_get(&timers_algos[B_KEY_GENERATION]), 2));
  std_timer[B_SIGN_ALGO] = sqrt(timer_pow2[B_SIGN_ALGO] - pow(btimer_get(&timers_algos[B_SIGN_ALGO]), 2));
  std_timer[B_VERIFY_ALGO] = sqrt(timer_pow2[B_VERIFY_ALGO] - pow(btimer_get(&timers_algos[B_VERIFY_ALGO]), 2));

  printf("Timing in ms:\n");
  printf(" - Key Gen: %.2f ms (std=%.2f)\n",
         btimer_get(&timers_algos[B_KEY_GENERATION]),
         std_timer[B_KEY_GENERATION]);
  printf(" - Sign:    %.2f ms (std=%.2f)\n",
         btimer_get(&timers_algos[B_SIGN_ALGO]),
         std_timer[B_SIGN_ALGO]);
  printf(" - Verify:  %.2f ms (std=%.2f)\n",
         btimer_get(&timers_algos[B_VERIFY_ALGO]),
         std_timer[B_VERIFY_ALGO]);
  printf("\n");

  printf("Timing in cycles:\n");
  printf(" - Key Gen: %.2f cycles\n", btimer_get_cycles(&timers_algos[B_KEY_GENERATION]));
  printf(" - Sign:    %.2f cycles\n", btimer_get_cycles(&timers_algos[B_SIGN_ALGO]));
  printf(" - Verify:  %.2f cycles\n", btimer_get_cycles(&timers_algos[B_VERIFY_ALGO]));
  printf("\n");

  free(workMem);
}
