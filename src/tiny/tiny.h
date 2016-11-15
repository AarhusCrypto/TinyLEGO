#ifndef TINY_TINY_TINY_H_
#define TINY_TINY_TINY_H_

#include "tiny/params.h"
#include "garbling/garbling-handler.h"
#include "circuit/circuit.h"

class Tiny {
public:
  Tiny(Params& params);

  void WeightedRandomString(uint8_t res[], int weight, int res_length, PRNG& rnd);

  virtual void Setup() = 0;
  virtual void Preprocess() = 0;
  virtual void Offline(std::vector<Circuit*>& circuits, int top_num_execs) = 0;
  
  Params& params;
  ctpl::thread_pool thread_pool;
  int num_gates_used;
  int num_inputs_used;
  int num_outputs_used;

  std::vector<int> gates_offset, inp_gates_offset, inputs_offset, outputs_offset;

  std::vector<std::unique_ptr<Params>> thread_params_vec;
  std::unique_ptr<uint8_t[]> thread_seeds;
};

#endif /* TINY_TINY_TINY_H_ */