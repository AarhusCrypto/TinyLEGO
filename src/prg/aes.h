// (C) 2016 University of Bristol. See LICENSE.txt

#ifndef __AES_H
#define __AES_H

#include "tiny-util/typedefs.h"

#define AES_BLK_SIZE 16

/************* C Version *************/
// Key Schedule
static inline void aes_schedule( uint* RK, uint8_t* K )
{ *RK = 0; *K = 0; }

// Encryption Function 
static inline void aes_encrypt( uint8_t* C, uint8_t* M, uint* RK )
{ *C = 0; *M = 0; *RK = 0; }


/*********** M-Code Version ***********/
// Check can support this
int Check_CPU_support_AES();
// Key Schedule 
void aes_128_schedule( uint8_t* key, const uint8_t* userkey );
void aes_192_schedule( uint8_t* key, const uint8_t* userkey );
void aes_256_schedule( uint8_t* key, const uint8_t* userkey );

static inline void aes_schedule( uint8_t* key, const uint8_t* userkey )
{ aes_128_schedule(key,userkey); }


// Encryption Function 
void aes_128_encrypt( uint8_t* C, const uint8_t* M,const uint8_t* RK );
void aes_192_encrypt( uint8_t* C, const uint8_t* M,const uint8_t* RK );
void aes_256_encrypt( uint8_t* C, const uint8_t* M,const uint8_t* RK );

static inline __m128i aes_128_encrypt(__m128i in, const uint8_t* key)
{ __m128i& tmp = in;
  tmp = _mm_xor_si128 (tmp,((__m128i*)key)[0]);
  int j;
  for(j=1; j <10; j++)
      { tmp = _mm_aesenc_si128 (tmp,((__m128i*)key)[j]); }
  tmp = _mm_aesenclast_si128 (tmp,((__m128i*)key)[j]);
  return tmp;
}

template <int N>

static inline void ecb_aes_128_encrypt(__m128i* out, __m128i* in, const uint8_t* key)
{
    __m128i tmp[N];
    for (int i = 0; i < N; i++)
        tmp[i] = _mm_xor_si128 (in[i],((__m128i*)key)[0]);
    int j;
    for(j=1; j <10; j++)
        for (int i = 0; i < N; i++)
            tmp[i] = _mm_aesenc_si128 (tmp[i],((__m128i*)key)[j]);
    for (int i = 0; i < N; i++)
        out[i] = _mm_aesenclast_si128 (tmp[i],((__m128i*)key)[j]);
}

template <int N>
static inline void ecb_aes_128_encrypt(__m128i* out, const __m128i* in, const uint8_t* key, const int* indices)
{
    __m128i tmp[N];
    for (int i = 0; i < N; i++)
        tmp[i] = in[indices[i]];
    ecb_aes_128_encrypt<N>(tmp, tmp, key);
    for (int i = 0; i < N; i++)
        out[indices[i]] = tmp[i];
}

static inline void aes_encrypt( uint8_t* C, const uint8_t* M,const uint8_t* RK )
{ aes_128_encrypt(C,M,RK); }

static inline __m128i aes_encrypt( __m128i M,const uint8_t* RK )
{ return aes_128_encrypt(M,RK); }


#endif

