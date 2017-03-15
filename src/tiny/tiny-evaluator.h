#ifndef TINY_TINY_TINYEVAL_H_
#define TINY_TINY_TINYEVAL_H_

#include "tiny/tiny.h"

// #include "dot/alsz-dot-ext-rec.h"
#include "commit/commit-scheme-rec.h"

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
  bool BatchDecommitLSB(CommitReceiver* commit_rec, uint8_t decommit_shares[], int num_values, uint8_t values[]);
  

  //For deletion soon
  // ALSZDOTExtRec ot_rec;
  // std::unique_ptr<uint8_t[]> rot_seeds;
  // std::unique_ptr<uint8_t[]> rot_choices;
  // std::vector<std::unique_ptr<CommitReceiver>> commit_recs;
  
  int rot_start_pos;

  std::unique_ptr<uint8_t[]> verleak_bits;
  std::unique_ptr<uint8_t[]> raw_eval_data;
  std::unique_ptr<uint32_t[]> raw_eval_ids;
  
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