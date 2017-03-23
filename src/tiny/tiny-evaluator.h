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
  void Online(std::vector<Circuit*>& circuits, std::vector<osuCrypto::BitVector>& inputs, std::vector<osuCrypto::BitVector>& outputs, int eval_num_execs);

  std::vector<uint8_t> raw_gates_data;
  std::vector<uint8_t> raw_auths_data;

  std::vector<uint32_t> eval_gates_ids;
  std::vector<uint32_t> eval_auths_ids;

  osuCrypto::BitVector global_dot_choices;
  osuCrypto::BitVector global_dot_lsb;
  BYTEArrayVector global_input_masks;

  BYTEArrayVector commit_shares_out_lsb_blind;
  
  std::vector<BYTEArrayVector> commit_shares;  

  std::vector<std::unique_ptr<osuCrypto::OtExtReceiver>> dot_receivers;
  std::vector<SplitCommitReceiver> commit_receivers;
  
  std::vector<osuCrypto::block> commit_seed_OTs;
  osuCrypto::BitVector commit_seed_choices;
  
  //Convenience pointers
  HalfGates eval_gates;
  Auths eval_auths;
};

#endif /* TINY_TINY_TINYEVAL_H_ */