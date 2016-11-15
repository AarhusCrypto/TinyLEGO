#include "test.h"

#include "commit/commit-scheme-snd.h"
#include "commit/commit-scheme-rec.h"

void RunSender(CommitSender& snd) {
  snd.Commit();
}

void RunReceiver(CommitReceiver& rec) {
  rec.Commit();
}

TEST(CommitCorrectness, Share0) {
  zmq::context_t context0(1);
  zmq::context_t context1(1);
  Params params_snd(constant_seeds[0], test_num_gates, test_num_inputs, test_num_outputs, default_ip_address, default_port, 0, context0, 2, GLOBAL_PARAMS_CHAN);
  Params params_rec(constant_seeds[1],  test_num_gates, test_num_inputs, test_num_outputs, default_ip_address, default_port, 1, context1, 2, GLOBAL_PARAMS_CHAN);

  std::unique_ptr<uint8_t[]> data(std::make_unique<uint8_t[]>(2 * CODEWORD_BITS * CSEC_BYTES + CODEWORD_BYTES));
  uint8_t* seeds0 = data.get();
  uint8_t* seeds1 = seeds0 + CODEWORD_BITS * CSEC_BYTES;
  uint8_t* choices = seeds1 + CODEWORD_BITS * CSEC_BYTES;
  
  PRNG rnd;
  rnd.SetSeed(constant_seeds[0]);
  rnd.GenRnd(seeds0, CODEWORD_BITS * CSEC_BYTES);
  rnd.GenRnd(seeds1, CODEWORD_BITS * CSEC_BYTES);

  CommitSender snd(params_snd, seeds0, seeds1);
  CommitReceiver rec(params_rec, seeds0, choices);

  mr_init_threading();
  thread snd_thread(RunSender, std::ref(snd));
  thread rec_thread(RunReceiver, std::ref(rec));
  snd_thread.join();
  rec_thread.join();
  mr_end_threading();

  for (int l = 0; l < params_snd.num_commits; l++) {
    for (int j = 0; j < CODEWORD_BITS; j++) {
      ASSERT_EQ(GetBitReversed(j, rec.commit_shares[l]), GetBitReversed(j, snd.commit_shares0[l]));
    }
  }
}

TEST(CommitCorrectness, Share1) {
  zmq::context_t context0(1);
  zmq::context_t context1(1);
  Params params_snd(constant_seeds[0], test_num_gates, test_num_inputs, test_num_outputs, default_ip_address, default_port, 0, context0, 2, GLOBAL_PARAMS_CHAN);
  Params params_rec(constant_seeds[1],  test_num_gates, test_num_inputs, test_num_outputs, default_ip_address, default_port, 1, context1, 2, GLOBAL_PARAMS_CHAN);

  std::unique_ptr<uint8_t[]> data(std::make_unique<uint8_t[]>(2 * CODEWORD_BITS * CSEC_BYTES + CODEWORD_BYTES));
  uint8_t* seeds0 = data.get();
  uint8_t* seeds1 = seeds0 + CODEWORD_BITS * CSEC_BYTES;
  uint8_t* choices = seeds1 + CODEWORD_BITS * CSEC_BYTES;
  std::fill(choices, choices + CODEWORD_BYTES, 0xFF);

  PRNG rnd;
  rnd.SetSeed(constant_seeds[0]);
  rnd.GenRnd(seeds0, CODEWORD_BITS * CSEC_BYTES);
  rnd.GenRnd(seeds1, CODEWORD_BITS * CSEC_BYTES);

  CommitSender snd(params_snd, seeds0, seeds1);
  CommitReceiver rec(params_rec, seeds1, choices);

  mr_init_threading();
  thread snd_thread(RunSender, std::ref(snd));
  thread rec_thread(RunReceiver, std::ref(rec));
  snd_thread.join();
  rec_thread.join();
  mr_end_threading();

  for (int l = 0; l < params_snd.num_commits; l++) {
    for (int j = 0; j < CODEWORD_BITS; j++) {
      ASSERT_EQ(GetBitReversed(j, rec.commit_shares[l]), GetBitReversed(j, snd.commit_shares1[l]));
    }
  }
}