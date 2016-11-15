#include "dot/alsz-dot-ext-rec.h"

ALSZDOTExtRec::ALSZDOTExtRec(Params& params) :
ALSZDOTExt(params),
choices_outer(std::make_unique<uint8_t[]>(BITS_TO_BYTES(params.num_OT))), 
response_outer(std::make_unique<uint8_t[]>(params.num_OT * CSEC_BYTES)),
receiver((crypto*)&params.crypt, net.rcvthread, net.sndthread, num_seed_OT, num_check_OT) {
}

void ALSZDOTExtRec::InitOTReceiver() {

  if (m_bUseMinEntCorAssumption) {
    receiver.EnableMinEntCorrRobustness();
  }
  receiver.ComputeBaseOTs(m_eFType);
}

void ALSZDOTExtRec::Receive() {
  CBitVector response_inner, choices;
  
  //Super hack. OTX wont work on too small inputs and how small depends on machine specs.
  int tmp_num_OT;
  if ((params.num_cpus == AWS_MACHINE_CORES) && (params.num_OT < AWS_MACHINE_MIN_OTX)) {
    tmp_num_OT = AWS_MACHINE_MIN_OTX; 
  } else if ((params.num_cpus == LLAN_MACHINE_CORES) && (params.num_OT < LLAN_MACHINE_MIN_OTX)) {
    tmp_num_OT = LLAN_MACHINE_MIN_OTX;
  } else {
    tmp_num_OT = params.num_OT;
  }
  
  response_inner.Create(tmp_num_OT, bit_length_inner);
  choices.Create(tmp_num_OT, (crypto*) &params.crypt);

  // Execute OT receiver routine
  auto OTX_begin = GET_TIME();
  receiver.receive(tmp_num_OT, bit_length_inner, num_snd_vals, &choices, &response_inner, s_type, r_type, num_OT_threads, m_fMaskFct.get());
  auto OTX_end = GET_TIME();
  
  //Then apply privacy amplification to response thus going from bit-strings of length k+2s to k. Can be k+s to k, but then we need so many checks in OTX that it takes longer than doing s more BaseOTs.
  auto privamp_begin = GET_TIME();
  PrivacyAmplification(response_inner.GetArr());
  std::copy(choices.GetArr(), choices.GetArr() + BITS_TO_BYTES(params.num_OT), choices_outer.get());
  auto privamp_end = GET_TIME();

  //Cleanup
  response_inner.delCBitVector();
  choices.delCBitVector();
#ifdef TINY_PRINT
  PRINT_TIME(OTX_end, OTX_begin, "OTX");
  PRINT_TIME(privamp_end, privamp_begin, "PRIVAMP");
#endif
}

void ALSZDOTExtRec::PrivacyAmplification(uint8_t response_inner[]) {
  //First receive the random seed from the constructor.
  std::unique_ptr<uint8_t[]> priv_amp_seed(std::make_unique<uint8_t[]>(CSEC_BYTES));
  params.chan.ReceiveBlocking(priv_amp_seed.get(), CSEC_BYTES);

  int byte_length_outer = BITS_TO_BYTES(bit_length_outer);
  uint8_t priv_amp_matrix[bit_length_inner * byte_length_outer];
  GeneratePrivAmpMatrix(priv_amp_seed.get(), priv_amp_matrix, bit_length_inner * byte_length_outer);

  ALSZDOTExt::PrivacyAmplification(priv_amp_matrix, byte_length_outer, bit_length_inner, params.num_OT, response_inner, response_outer.get());
}