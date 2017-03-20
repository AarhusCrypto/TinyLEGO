#include "garbling/garbling-handler.h"

GarblingHandler::GarblingHandler(Params& params) : params(params) {
  aes128_load_key(global_aes_key, key_schedule);
}
//Below commented out part can be used for debugging in tiny-constructor and tiny-evaluator
/////////////////// DEBUG for testing correctness of solderings//////////////
bool GarblingHandler::AllVerifyAuths(Auths& auths_data, int offset, uint8_t keys[], uint32_t ids[], uint32_t num_auths) {
  __m128i key_128, hash_128, id_128, S_A_128;
  uint64_t tmp_id;
  for (uint32_t i = 0; i < num_auths; ++i) {
    int current_auth_index = offset + i;

    key_128 = _mm_lddqu_si128((__m128i *) (keys)); 

    hash_128 = _mm_lddqu_si128((__m128i *) (auths_data.H_0 + current_auth_index * CSEC_BYTES));
    tmp_id = ids[current_auth_index]; //Upcast to 64bit for loadl
    id_128 = _mm_loadl_epi64((__m128i *) &tmp_id);

    //Soldering
    S_A_128 = _mm_lddqu_si128((__m128i *) (auths_data.S_A + current_auth_index * CSEC_BYTES));

    key_128 = _mm_xor_si128(key_128, S_A_128);

    key_128 = AESHash(key_128, id_128);
    if (compare128(key_128, hash_128)) {
      //Matched first authenticator
    } else {
      hash_128 = _mm_lddqu_si128((__m128i *) (auths_data.H_1 + current_auth_index * CSEC_BYTES));
      if (compare128(key_128, hash_128)) {
        //Matched second authenticator
      } else {
        //Matched none of the authenticators
        return false;
      }
    }
  }
  return true;
}

void GarblingHandler::AllShiftEvaluateGates(HalfGates& gates_data, int offset, uint8_t left_keys[], uint8_t right_keys[], uint8_t out_keys[], uint32_t ids[], uint32_t num_gates) {
  __m128i left_key_128, right_key_128, id_128, T_G_128, T_E_128, out_key_128, tmp_128, S_L_128, S_R_128, S_O_128;
  uint64_t tmp_id;
  for (uint32_t i = 0; i < num_gates; ++i) {
    int current_gate_index = offset + i;

    left_key_128 = _mm_lddqu_si128((__m128i *) (left_keys));
    right_key_128 = _mm_lddqu_si128((__m128i *) (right_keys));
    T_G_128 = _mm_lddqu_si128((__m128i *) (gates_data.T_G + current_gate_index * CSEC_BYTES));
    T_E_128 = _mm_lddqu_si128((__m128i *) (gates_data.T_E + current_gate_index * CSEC_BYTES));
    tmp_id = ids[current_gate_index]; //Upcast to 64bit for loadl
    id_128 = _mm_loadl_epi64((__m128i *) &tmp_id);

    //Soldering
    S_L_128 = _mm_lddqu_si128((__m128i *) (gates_data.S_L + current_gate_index * CSEC_BYTES));
    S_R_128 = _mm_lddqu_si128((__m128i *) (gates_data.S_R + current_gate_index * CSEC_BYTES));
    S_O_128 = _mm_lddqu_si128((__m128i *) (gates_data.S_O + current_gate_index * CSEC_BYTES));


    left_key_128 = _mm_xor_si128(left_key_128, S_L_128);
    right_key_128 = _mm_xor_si128(right_key_128, S_R_128);

    out_key_128 = AESHash(left_key_128, id_128);
    tmp_128 = AESHash(right_key_128, id_128);

    if (GetLSB(left_key_128)) {
      out_key_128 = _mm_xor_si128(out_key_128, T_G_128);
    }

    out_key_128 = _mm_xor_si128(out_key_128, tmp_128);

    if (GetLSB(right_key_128)) {
      out_key_128 = _mm_xor_si128(out_key_128, T_E_128);
      out_key_128 = _mm_xor_si128(out_key_128, left_key_128);
    }

    out_key_128 = _mm_xor_si128(out_key_128, S_O_128);
    _mm_storeu_si128((__m128i *) (out_keys), out_key_128);
  }
}

//HalfGate garbling
void GarblingHandler::GarbleGates(HalfGates& gates_data, int offset, uint8_t left_keys[], uint8_t right_keys[], uint8_t delta[], uint32_t ids[], uint32_t num_gates) {
  __m128i delta_128 = _mm_lddqu_si128((__m128i *) delta);
  __m128i left_key_128, right_key_128, left_key_delta_128, right_key_delta_128, out_key_128, T_G_128, T_E_128, id_128, tmp_128;
  for (uint32_t i = 0; i < num_gates; ++i) {
    uint8_t left_bit = GetLSB(left_keys + (offset + i) * CSEC_BYTES);
    uint8_t right_bit = GetLSB(right_keys + (offset + i) * CSEC_BYTES);

    left_key_128 = _mm_lddqu_si128((__m128i *) (left_keys + (offset + i) * CSEC_BYTES));
    right_key_128 = _mm_lddqu_si128((__m128i *) (right_keys + (offset + i) * CSEC_BYTES));
    id_128 = (__m128i) _mm_load_ss((float*) &ids[offset + i]);

    left_key_delta_128 = _mm_xor_si128(left_key_128, delta_128);
    right_key_delta_128 = _mm_xor_si128(right_key_128, delta_128);

    out_key_128 = AESHash(left_key_128, id_128);
    T_G_128 = AESHash(left_key_delta_128, id_128);
    T_E_128 = AESHash(right_key_128, id_128);
    tmp_128 = AESHash(right_key_delta_128, id_128);

    T_G_128 = _mm_xor_si128(T_G_128, out_key_128);
    if (right_bit) {
      T_G_128 = _mm_xor_si128(T_G_128, delta_128);
      out_key_128 = _mm_xor_si128(out_key_128, tmp_128);
    } else {
      out_key_128 = _mm_xor_si128(out_key_128, T_E_128);
    }

    if (left_bit) {
      out_key_128 = _mm_xor_si128(out_key_128, T_G_128);
    }

    T_E_128 = _mm_xor_si128(T_E_128, left_key_128);
    T_E_128 = _mm_xor_si128(T_E_128, tmp_128);

    _mm_storeu_si128((__m128i *) (gates_data.T_G + (offset + i) * CSEC_BYTES), T_G_128);
    _mm_storeu_si128((__m128i *) (gates_data.T_E + (offset + i) * CSEC_BYTES), T_E_128);
    _mm_storeu_si128((__m128i *) (gates_data.S_O + (offset + i) * CSEC_BYTES), out_key_128);
  }
}

//HalfGate Evaluation
void GarblingHandler::OutputShiftEvaluateGates(HalfGates& gates_data, int offset, uint8_t left_keys[], uint8_t right_keys[], uint8_t out_keys[], uint32_t ids[], uint32_t num_gates, int neg_offset_ids) {
  __m128i left_key_128, right_key_128, id_128, T_G_128, T_E_128, out_key_128, tmp_128, S_L_128, S_R_128, S_O_128;
  for (uint32_t i = 0; i < num_gates; ++i) {
    uint8_t left_bit = GetLSB(left_keys + (offset + i) * CSEC_BYTES);
    uint8_t right_bit = GetLSB(right_keys + (offset + i) * CSEC_BYTES);

    int current_gate_index = ids[offset + i] - neg_offset_ids - params.out_keys_start;

    S_O_128 = _mm_lddqu_si128((__m128i *) (gates_data.S_O + current_gate_index * CSEC_BYTES));

    left_key_128 = _mm_lddqu_si128((__m128i *) (left_keys + (offset + i) * CSEC_BYTES));
    right_key_128 = _mm_lddqu_si128((__m128i *) (right_keys + (offset + i) * CSEC_BYTES));
    id_128 = (__m128i) _mm_load_ss((float*) &ids[offset + i]);
    T_G_128 = _mm_lddqu_si128((__m128i *) (gates_data.T_G + current_gate_index * CSEC_BYTES));
    T_E_128 = _mm_lddqu_si128((__m128i *) (gates_data.T_E + current_gate_index * CSEC_BYTES));

    out_key_128 = AESHash(left_key_128, id_128);
    tmp_128 = AESHash(right_key_128, id_128);

    if (left_bit) {
      out_key_128 = _mm_xor_si128(out_key_128, T_G_128);
    }

    out_key_128 = _mm_xor_si128(out_key_128, tmp_128);
    if (right_bit) {
      out_key_128 = _mm_xor_si128(out_key_128, T_E_128);
      out_key_128 = _mm_xor_si128(out_key_128, left_key_128);
    }

    out_key_128 = _mm_xor_si128(out_key_128, S_O_128);
    _mm_storeu_si128((__m128i *) (out_keys + (offset + i) * CSEC_BYTES), out_key_128);
  }
}

//Wire Authenticators production
void GarblingHandler::GarbleAuths(Auths& auths_data, int offset, uint8_t keys[], uint8_t delta[], uint32_t ids[], uint32_t num_auths) {
  __m128i delta_128 = _mm_lddqu_si128((__m128i *) delta);
  __m128i key_128, key_delta_128, id_128;
  for (uint32_t i = 0; i < num_auths; ++i) {
    key_128 = _mm_lddqu_si128((__m128i *) (keys + (offset + i) * CSEC_BYTES));
    id_128 = (__m128i) _mm_load_ss((float*) &ids[offset + i]);

    key_delta_128 = _mm_xor_si128(key_128, delta_128);
    key_128 = AESHash(key_128, id_128);
    key_delta_128 = AESHash(key_delta_128, id_128);
    _mm_storeu_si128((__m128i *) (auths_data.H_0 + (offset + i) * CSEC_BYTES), key_128);

    _mm_storeu_si128((__m128i *) (auths_data.H_1 + (offset + i) * CSEC_BYTES), key_delta_128);
    int res = memcmp(auths_data.H_0 + (offset + i) * CSEC_BYTES, auths_data.H_1 + (offset + i) * CSEC_BYTES, CSEC_BYTES);
    if (res > 0) {
      //Do nothing
    } else if (res < 0) {
      //Swap the order of the authenticators
      uint8_t tmp[CSEC_BYTES];
      std::copy(auths_data.H_0 + (offset + i) * CSEC_BYTES, auths_data.H_0 + (offset + i) * CSEC_BYTES + CSEC_BYTES, tmp);
      std::copy(auths_data.H_1 + (offset + i) * CSEC_BYTES, auths_data.H_1 + (offset + i) * CSEC_BYTES + CSEC_BYTES, auths_data.H_0 + (offset + i) * CSEC_BYTES);
      std::copy(tmp, tmp + CSEC_BYTES, auths_data.H_1 + (offset + i) * CSEC_BYTES);
    } else {
      std::cout << "Congrats, this only happens with prob. 2^-128! It must be your lucky day!" << std::endl;
    }
  }
}

//Verify Wire Authenticators
bool GarblingHandler::VerifyAuths(Auths& auths_data, int offset, uint8_t keys[], uint32_t ids[], uint32_t num_auths, int neg_offset_ids) {
  __m128i key_128, hash_128, id_128;
  for (uint32_t i = 0; i < num_auths; ++i) {
    int current_auth_index = ids[offset + i] - neg_offset_ids - params.auth_start;

    key_128 = _mm_lddqu_si128((__m128i *) (keys + (offset + i) * CSEC_BYTES));
    hash_128 = _mm_lddqu_si128((__m128i *) (auths_data.H_0 + current_auth_index * CSEC_BYTES));
    id_128 = (__m128i) _mm_load_ss((float*) &ids[offset + i]);
    key_128 = AESHash(key_128, id_128);
    if (compare128(key_128, hash_128)) {
      //Matched first authenticator
    } else {
      hash_128 = _mm_lddqu_si128((__m128i *) (auths_data.H_1 + current_auth_index * CSEC_BYTES));
      if (compare128(key_128, hash_128)) {
        //Matched second authenticator
      } else {
        //Matched none of the authenticators
        return false;
      }
    }
  }
  return true;
}

//Fixed-Key AES Hash
__m128i GarblingHandler::AESHash(__m128i& value_128, __m128i& id_128) {
  __m128i res = DOUBLE(value_128);
  res = _mm_xor_si128(res, id_128);
  __m128i tmp_input = res;

  DO_ENC_BLOCK(res, key_schedule);

  return _mm_xor_si128(res, tmp_input);
}