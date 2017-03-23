#ifndef TINY_TINY_TINYCONST_H_
#define TINY_TINY_TINYCONST_H_

#include "tiny/tiny.h"

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
  void Online(std::vector<Circuit*>& circuits, std::vector<osuCrypto::BitVector>& inputs, int eval_num_execs);

  std::vector<uint32_t> eval_gates_ids;
  std::vector<uint32_t> eval_auths_ids;

  std::vector<std::array<BYTEArrayVector, 2>> commit_shares;
  std::vector<std::array<osuCrypto::block, 2>> commit_seed_OTs;
  std::vector<std::unique_ptr<osuCrypto::OtExtSender>> dot_senders;
  std::vector<SplitCommitSender> commit_senders;

  std::array<BYTEArrayVector, 2> commit_shares_out_lsb_blind;

};

#endif /* TINY_TINY_TINYCONST_H_ */