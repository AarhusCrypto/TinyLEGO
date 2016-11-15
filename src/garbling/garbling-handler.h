#ifndef TINY_GARBLING_GARBLINGHANDLER_H_
#define TINY_GARBLING_GARBLINGHANDLER_H_

#include "garbling/halfgates.h"
#include "garbling/auths.h"

#include "tiny/params.h"

class GarblingHandler {
public:
  GarblingHandler(Params& params);

  void AllShiftEvaluateGates(HalfGates& gates_data, int offset, uint8_t left_keys[], uint8_t right_keys[], uint8_t out_keys[], uint32_t ids[], uint32_t num_gates);
  bool AllVerifyAuths(Auths& auths_data, int offset, uint8_t keys[], uint32_t ids[], uint32_t num_auths);

  void GarbleGates(HalfGates& gates_data, int offset, uint8_t left_keys[], uint8_t right_keys[], uint8_t delta[], uint32_t ids[], uint32_t num_gates);

  void OutputShiftEvaluateGates(HalfGates& gates_data, int offset, uint8_t left_keys[], uint8_t right_keys[], uint8_t out_keys[], uint32_t ids[], uint32_t num_gates, int neg_offset_ids);

  void GarbleAuths(Auths& auths_data, int offset, uint8_t keys[], uint8_t delta[], uint32_t ids[], uint32_t num_auths);

  bool VerifyAuths(Auths& auths_data, int offset, uint8_t keys[], uint32_t ids[], uint32_t num_auths, int neg_offset_ids);

  __m128i AESHash(__m128i& value_128, __m128i& id_128);

  Params& params;
  __m128i key_schedule[11];
};

//Static context for added performance
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

static inline __m128i aes_128_key_expansion(__m128i key, __m128i keygened) {
  keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3, 3, 3, 3));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  return _mm_xor_si128(key, keygened);
}

#define AES_128_key_exp(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))

static inline void aes128_load_key(uint8_t *enc_key, __m128i key_schedule[]) {
  key_schedule[0]  = _mm_lddqu_si128((const __m128i*) enc_key);
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
};

//HalfGate Evaluation
static inline void IntrinShiftEvaluateGates(HalfGates& gates_data, int offset, __m128i& left_key_128, __m128i& right_key_128, __m128i& out_key_128, uint32_t id, __m128i key_schedule[]) {

  __m128i T_G_128 = _mm_lddqu_si128((__m128i *) (gates_data.T_G + offset * AES_BYTES));
  __m128i T_E_128 = _mm_lddqu_si128((__m128i *) (gates_data.T_E + offset * AES_BYTES));
  __m128i id_128 = (__m128i) _mm_load_ss((float*) &id);

  //Soldering
  __m128i S_L_128 = _mm_lddqu_si128((__m128i *) (gates_data.S_L + offset * AES_BYTES));
  __m128i S_R_128 = _mm_lddqu_si128((__m128i *) (gates_data.S_R + offset * AES_BYTES));
  __m128i S_O_128 = _mm_lddqu_si128((__m128i *) (gates_data.S_O + offset * AES_BYTES));
  S_L_128 = _mm_xor_si128(left_key_128, S_L_128);
  S_R_128 = _mm_xor_si128(right_key_128, S_R_128);

  //////First Hash//////
  __m128i res = _mm_xor_si128(DOUBLE(S_L_128), id_128);
  // __m128i out_key_128 = res;
  out_key_128 = res;
  DO_ENC_BLOCK(res, key_schedule);
  out_key_128 = _mm_xor_si128(out_key_128, res);
  //////First Hash//////

  //////Second Hash//////
  res = _mm_xor_si128(DOUBLE(S_R_128), id_128);
  out_key_128 = _mm_xor_si128(out_key_128, res);
  DO_ENC_BLOCK(res, key_schedule);
  out_key_128 = _mm_xor_si128(out_key_128, res);
  //////Second Hash//////

  out_key_128 = _mm_xor_si128(out_key_128, _mm_and_si128(T_G_128, invert_array[GetLSB(S_L_128)]));
  //Equals to
  // if (GetLSB(S_L_128)) {
  // out_key_128 = _mm_xor_si128(out_key_128, T_G_128);
  // }

  out_key_128 = _mm_xor_si128(out_key_128, _mm_and_si128(T_E_128, invert_array[GetLSB(S_R_128)]));
  out_key_128 = _mm_xor_si128(out_key_128, _mm_and_si128(S_L_128, invert_array[GetLSB(S_R_128)]));
  //Equals to
  // if (GetLSB(S_R_128)) {
  //   out_key_128 = _mm_xor_si128(out_key_128, T_E_128);
  //   out_key_128 = _mm_xor_si128(out_key_128, S_L_128);
  // }

  //Output soldering
  out_key_128 = _mm_xor_si128(out_key_128, S_O_128);

};

static inline bool IntrinVerifyAuths(Auths& auths_data, int offset, __m128i key_128, uint32_t id, __m128i key_schedule[]) {
  // __m128i hash_128, id_128, S_A_128;

  __m128i hash_128 = _mm_lddqu_si128((__m128i *) (auths_data.H_0 + offset * AES_BYTES));
  __m128i id_128 = (__m128i) _mm_load_ss((float*) &id);

  //Soldering
  __m128i S_A_128 = _mm_lddqu_si128((__m128i *) (auths_data.S_A + offset * AES_BYTES));

  key_128 = _mm_xor_si128(key_128, S_A_128);

  // key_128 = AESHash(key_128, id_128);
  __m128i res = _mm_xor_si128(DOUBLE(key_128), id_128);
  key_128 = res;
  DO_ENC_BLOCK(res, key_schedule);
  key_128 = _mm_xor_si128(key_128, res);

  if (compare128(key_128, hash_128)) {
    //Matched first authenticator
  } else {
    hash_128 = _mm_lddqu_si128((__m128i *) (auths_data.H_1 + offset * AES_BYTES));
    if (compare128(key_128, hash_128)) {
      //Matched second authenticator
    } else {
      //Matched none of the authenticators
      return false;
    }
  }

  return true;
};

#endif /* TINY_GARBLING_GARBLINGHANDLER_H_ */