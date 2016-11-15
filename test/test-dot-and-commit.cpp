#include "test.h"

#include "commit/commit-scheme-snd.h"
#include "commit/commit-scheme-rec.h"
#include "dot/alsz-dot-ext-snd.h"
#include "dot/alsz-dot-ext-rec.h"
#include "tiny/params.h"

void RunSenderDOT(ALSZDOTExtSnd& dot_snd) {
  dot_snd.InitOTSender();
  dot_snd.Send();
}

void RunReceiverDOT(ALSZDOTExtRec& dot_rec) {
  dot_rec.InitOTReceiver();
  dot_rec.Receive();
}

void RunSenderCommit(CommitSender& commit_snd) {
  commit_snd.Commit();
}

void RunReceiverCommit(CommitReceiver& commit_rec) {
  commit_rec.Commit();
}

TEST(DOTAndCommit, Full) {
  zmq::context_t context0(1);
  zmq::context_t context1(1);
  Params params_snd(constant_seeds[0], test_num_gates, test_num_inputs, test_num_outputs, default_ip_address, default_port, 0, context0, 2, GLOBAL_PARAMS_CHAN);
  Params params_rec(constant_seeds[1],  test_num_gates, test_num_inputs, test_num_outputs, default_ip_address, default_port, 1, context1, 2, GLOBAL_PARAMS_CHAN);

  ALSZDOTExtSnd dot_snd(params_snd, true);
  ALSZDOTExtRec dot_rec(params_rec);

  mr_init_threading();
  thread snd_dot_thread(RunSenderDOT, std::ref(dot_snd));
  thread rec_dot_thread(RunReceiverDOT, std::ref(dot_rec));

  snd_dot_thread.join();
  rec_dot_thread.join();
  mr_end_threading();

  //Construct base[i] xor Delta
  std::unique_ptr<uint8_t[]> base_delta_outer_ptr(std::make_unique<uint8_t[]>(params_snd.num_OT * CSEC_BYTES));
  uint8_t* base_delta_outer = base_delta_outer_ptr.get();

  for (int i = 0; i < params_snd.num_OT; i++) {
    XOR_128(base_delta_outer + i * CSEC_BYTES, dot_snd.base_outer.get() + i * CSEC_BYTES, dot_snd.delta_outer.get());
  }

  for (int j = 0; j < CODEWORD_BITS; j++) {
    if (GetBit(j, dot_rec.choices_outer.get())) {
      for (int i = 0; i < CSEC_BYTES; i++) {
        ASSERT_EQ(base_delta_outer[j * CSEC_BYTES + i], dot_rec.response_outer.get()[j * CSEC_BYTES + i]);
      }
    } else {
      for (int i = 0; i < CSEC_BYTES; i++) {
        ASSERT_EQ(dot_snd.base_outer.get()[j * CSEC_BYTES + i], dot_rec.response_outer.get()[j * CSEC_BYTES + i]);
      }
    }
  }

  CommitSender commit_snd(params_snd, dot_snd.base_outer.get(), base_delta_outer);
  CommitReceiver commit_rec(params_rec, dot_rec.response_outer.get(), dot_rec.choices_outer.get());

  mr_init_threading();
  thread snd_commit_thread(RunSenderCommit, std::ref(commit_snd));
  thread rec_commit_thread(RunReceiverCommit, std::ref(commit_rec));
  snd_commit_thread.join();
  rec_commit_thread.join();
  mr_end_threading();

  for (int l = 0; l < params_snd.num_commits; l++) {
    for (int j = 0; j < CODEWORD_BITS; j++) {
      if (GetBit(j, dot_rec.choices_outer.get())) {
        ASSERT_EQ(GetBitReversed(j, commit_rec.commit_shares[l]), GetBitReversed(j, commit_snd.commit_shares1[l]));
      } else {
        ASSERT_EQ(GetBitReversed(j, commit_rec.commit_shares[l]), GetBitReversed(j, commit_snd.commit_shares0[l]));
      }
    }
  }
}