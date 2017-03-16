#include "tiny/tiny.h"

Tiny::Tiny(uint8_t seed[], Params& params) :
params(params),
thread_pool(params.num_cpus),
exec_rnds(params.num_execs),
ios(0) {
  
  rnd.SetSeed(load_block(seed));

  uint8_t tmp_seed[CSEC_BYTES];
  for (int e = 0; e < params.num_execs; ++e) {
    rnd.get<uint8_t>(tmp_seed, CSEC_BYTES);
    exec_rnds[e].SetSeed(load_block(tmp_seed));
  }

  //Makes bucketing much easier to handle. The AES_BITS factor is needed to order to make VerLeak check behave properly
  int num_gates_rounded = PAD_TO_MULTIPLE(params.num_pre_gates, params.num_execs);
  int num_inputs_rounded = PAD_TO_MULTIPLE(params.num_pre_inputs, AES_BITS * params.num_execs);
  int num_outputs_rounded = PAD_TO_MULTIPLE(params.num_pre_outputs, AES_BITS * params.num_execs);

  params.ComputeGateAndAuthNumbers(num_gates_rounded, num_inputs_rounded, num_outputs_rounded);

  num_gates_used = 0;
  num_inputs_used = 0;
  num_outputs_used = 0;
}

Tiny::~Tiny() {
  ios.stop();
}

void Tiny::WeightedRandomString(uint8_t res[], int weight, int res_length, osuCrypto::PRNG& rnd) {

  uint8_t temp[res_length];
  std::fill(res, res + res_length, 0xFF);
  while (weight > 0) {
    rnd.get<uint8_t>(temp, res_length);
    for (unsigned int i = 0; i < res_length; ++i) {
      res[i] &= temp[i];
    }
    weight--;
  }
}