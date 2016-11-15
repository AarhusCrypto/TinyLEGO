#include "dot/alsz-dot-ext-snd.h"

ALSZDOTExtSnd::ALSZDOTExtSnd(Params& params, bool set_lsb_delta) :
  ALSZDOTExt(params),
  set_lsb_delta(set_lsb_delta),
  base_outer(std::make_unique<uint8_t[]>(params.num_OT * CSEC_BYTES)),
  delta_outer(std::make_unique<uint8_t[]>(CSEC_BYTES)),
  sender(ALSZOTExtSnd((crypto*) & params.crypt, net.rcvthread, net.sndthread, num_seed_OT, num_check_OT)) {
}

void ALSZDOTExtSnd::InitOTSender() {

  if (m_bUseMinEntCorAssumption) {
    sender.EnableMinEntCorrRobustness();
  }
  sender.ComputeBaseOTs(m_eFType);
}

void ALSZDOTExtSnd::Send() {

  int byte_length_inner = BITS_TO_BYTES(bit_length_inner);
  uint8_t delta_inner[byte_length_inner];

  //Super hack. OTX wont work on too small inputs and how small depends on machine specs.
  int tmp_num_OT;
  if ((params.num_cpus == AWS_MACHINE_CORES) && (params.num_OT < AWS_MACHINE_MIN_OTX)) {
    tmp_num_OT = AWS_MACHINE_MIN_OTX; 
  } else if ((params.num_cpus == LLAN_MACHINE_CORES) && (params.num_OT < LLAN_MACHINE_MIN_OTX)) {
    tmp_num_OT = LLAN_MACHINE_MIN_OTX;
  } else {
    tmp_num_OT = params.num_OT;
  }

  //Cannot use unique_ptr due to interface of OTX
  CBitVector** X = new CBitVector*[num_snd_vals]; //Ownership is passed on to sender, so it gets deleted when sender is deleted.
  X[0] = new CBitVector();
  X[1] = new CBitVector();

  X[0]->Create(tmp_num_OT * bit_length_inner);
  X[1]->Create(tmp_num_OT * bit_length_inner);

  // Execute OT sender routine
  auto OTX_begin = GET_TIME();
  sender.send(tmp_num_OT, bit_length_inner, num_snd_vals, X, s_type, r_type, num_OT_threads, m_fMaskFct.get());
  auto OTX_end = GET_TIME();

  //Post OTX processing
  std::copy(X[1]->GetArr(), X[1]->GetArr() + byte_length_inner, delta_inner);
  XOR_UINT8_T(delta_inner, X[0]->GetArr(), byte_length_inner);


  //Then apply privacy amplification to vectors and delta thus going from bit-strings of length k+s to k. Can be k+s to k, but then we need so many checks in OTX that it takes longer than doing s more BaseOTs.
  auto privamp_begin = GET_TIME();
  PrivacyAmplification(X[0]->GetArr(), delta_inner);
  auto privamp_end = GET_TIME();

  X[0]->delCBitVector();
  X[1]->delCBitVector();
  delete X[0];
  delete X[1];

#ifdef TINY_PRINT
  PRINT_TIME(OTX_end, OTX_begin, "OTX");
  PRINT_TIME(privamp_end, privamp_begin, "PRIVAMP");
#endif
}

void ALSZDOTExtSnd::PrivacyAmplification(uint8_t base_inner[], uint8_t delta_inner[]) {

  int byte_length_outer = BITS_TO_BYTES(bit_length_outer);
  uint8_t priv_amp_seed[CSEC_BYTES];
  uint8_t priv_amp_matrix[bit_length_inner * byte_length_outer];

  bool done = false;
  while (!done) {
    params.crypt.gen_rnd(priv_amp_seed, CSEC_BYTES);
    GeneratePrivAmpMatrix(priv_amp_seed, priv_amp_matrix, bit_length_inner * byte_length_outer);
    for (int bit = 0; bit < bit_length_inner; ++bit) {
      if (GetBitReversed(bit, delta_inner)) {
        XOR_128(delta_outer.get(), priv_amp_matrix + (bit * byte_length_outer));
      }
    }
    //If set_lsb_delta flag is set, this ensures that lsb(delta) == 1. This is needed for Half-Gate garbling.
    if (set_lsb_delta && (GetLSB(delta_outer.get()) != 1)) {
      //Reset delta_out as we're going to go into the loop again
      std::fill(delta_outer.get(), delta_outer.get() + CSEC_BYTES, 0);
    } else {
      //We exit loop
      done = true;
    }
  }

  params.chan.Send(priv_amp_seed, CSEC_BYTES);

  ALSZDOTExt::PrivacyAmplification(priv_amp_matrix, byte_length_outer, bit_length_inner, params.num_OT, base_inner, base_outer.get());
}