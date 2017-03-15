#ifndef TINY_TINY_TINY_H_
#define TINY_TINY_TINY_H_

#include "tiny/params.h"
#include "garbling/garbling-handler.h"
#include "circuit/circuit.h"

//For libOTe extension
#include "cryptoTools/Network/BtChannel.h"
#include "cryptoTools/Network/BtEndpoint.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/BitVector.h"
#include "libOTe/Base/naor-pinkas.h"

#include "cryptoTools/Common/Defines.h"

class Tiny {
public:
  Tiny(uint8_t seed[], Params& params);
  ~Tiny();

  void WeightedRandomString(uint8_t res[], int weight, int res_length, osuCrypto::PRNG& rnd);

  virtual void Setup() = 0;
  virtual void Preprocess() = 0;
  virtual void Offline(std::vector<Circuit*>& circuits, int top_num_execs) = 0;
  
  Params& params;
  ctpl::thread_pool thread_pool;
  int num_gates_used;
  int num_inputs_used;
  int num_outputs_used;

  std::vector<int> gates_offset, inp_gates_offset, inputs_offset, outputs_offset;

  osuCrypto::PRNG rnd;

  std::vector<osuCrypto::PRNG> exec_rnds;
  osuCrypto::BtIOService ios;
  osuCrypto::BtEndpoint end_point;
  osuCrypto::Channel* chan;
  std::vector<osuCrypto::Channel*> exec_channels;

  std::vector<std::unique_ptr<Params>> thread_params_vec;
  std::unique_ptr<uint8_t[]> thread_seeds;
};

#endif /* TINY_TINY_TINY_H_ */