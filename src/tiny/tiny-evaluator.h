#ifndef TINY_TINY_TINYEVAL_H_
#define TINY_TINY_TINYEVAL_H_

#include "tiny/tiny.h" 

#include "libOTe/TwoChooseOne/KosOtExtReceiver.h"
#include "libOTe/TwoChooseOne/KosDotExtReceiver.h"
#include "split-commit/split-commit-rec.h"

class TinyEvaluator : public Tiny {
public:
  TinyEvaluator(uint8_t seed[], Params& params);

  ~TinyEvaluator();
  void Connect(std::string ip_address, uint16_t port);

  void Setup();
  void Preprocess();
  void Offline(std::vector<Circuit*>& circuits, int top_num_execs);
  void Online(std::vector<Circuit*>& circuits, std::vector<uint8_t*>& inputs, std::vector<uint8_t*>& outputs, int eval_num_execs);

  std::unique_ptr<uint8_t[]> raw_eval_data;
  std::unique_ptr<uint32_t[]> raw_eval_ids;

  osuCrypto::BitVector global_dot_choices;
  osuCrypto::BitVector global_dot_lsb;
  osuCrypto::BitVector global_out_lsb;
  BYTEArrayVector global_input_masks;
  
  std::vector<BYTEArrayVector> commit_shares;  

  std::vector<std::unique_ptr<osuCrypto::OtExtReceiver>> dot_receivers;
  std::vector<SplitCommitReceiver> commit_receivers;
  
  std::vector<osuCrypto::block> commit_seed_OTs;
  osuCrypto::BitVector commit_seed_choices;
  
  //Convenience pointers
  uint32_t* eval_gates_ids;
  uint32_t* eval_auths_ids;
  HalfGates eval_gates;
  Auths eval_auths;
};

#endif /* TINY_TINY_TINYEVAL_H_ */