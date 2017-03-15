#ifndef TINY_TINY_TINYCONST_H_
#define TINY_TINY_TINYCONST_H_

#include "tiny/tiny.h"

// #include "dot/alsz-dot-ext-snd.h"
#include "commit/commit-scheme-snd.h"

#include "libOTe/TwoChooseOne/KosOtExtSender.h"
#include "libOTe/TwoChooseOne/KosDotExtSender.h"
#include "split-commit/split-commit-snd.h"

class TinyConstructor : public Tiny {
public:
  TinyConstructor(uint8_t seed[], Params& params);

  ~TinyConstructor();

  void Connect(std::string ip_address, uint16_t port);

  void Setup();
  void Preprocess();
  void Offline(std::vector<Circuit*>& circuits, int top_num_execs);
  void Online(std::vector<Circuit*>& circuits, std::vector<uint8_t*>& inputs, int eval_num_execs);
  void BatchDecommitLSB(CommitSender* commit_snd, uint8_t decommit_shares0[], uint8_t decommit_shares1[], int num_values);

  //For deletion soon
  // ALSZDOTExtSnd ot_snd;
  // std::unique_ptr<uint8_t[]> rot_seeds0;
  // uint8_t* rot_seeds1;
  // std::vector<std::unique_ptr<CommitSender>> commit_snds;

  std::unique_ptr<uint32_t[]> raw_eval_ids;

  std::vector<std::array<BYTEArrayVector, 2>> commit_shares;
  std::vector<std::array<osuCrypto::block, 2>> commit_seed_OTs;
  std::vector<std::unique_ptr<osuCrypto::OtExtSender>> dot_senders;
  std::vector<SplitCommitSender> commit_senders;


  uint32_t* eval_gates_ids;
  uint32_t* eval_auths_ids;
  uint8_t global_delta[CSEC_BYTES] = {0};
};

#endif /* TINY_TINY_TINYCONST_H_ */