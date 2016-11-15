#include "test.h"

#include "dot/alsz-dot-ext-snd.h"
#include "dot/alsz-dot-ext-rec.h"

void RunSender(ALSZDOTExtSnd& snd) {
  snd.InitOTSender();
  snd.Send();
  
}

void RunReceiver(ALSZDOTExtRec& rec) {
  rec.InitOTReceiver();
  rec.Receive();
}


TEST(FULL_DOT, Correctness) {
  zmq::context_t context0(1);
  zmq::context_t context1(1);
  Params params_snd(constant_seeds[0], test_num_gates, test_num_inputs, test_num_outputs, default_ip_address, default_port, 0, context0, 2, GLOBAL_PARAMS_CHAN);
  Params params_rec(constant_seeds[1],  test_num_gates, test_num_inputs, test_num_outputs, default_ip_address, default_port, 1, context1, 2, GLOBAL_PARAMS_CHAN);
  
  ALSZDOTExtSnd snd(params_snd, true);
  ALSZDOTExtRec rec(params_rec);

  mr_init_threading();
  thread snd_thread(RunSender, std::ref(snd));
  thread rec_thread(RunReceiver, std::ref(rec));
  snd_thread.join();
  rec_thread.join();
  mr_end_threading();

  uint8_t tmp[CSEC_BYTES] = {0};
  for (int i = 0; i < params_snd.num_OT; i++) {
    if (GetBit(i, rec.choices_outer.get())) {
      XOR_128(tmp, snd.base_outer.get() + i * CSEC_BYTES, snd.delta_outer.get());
      for (int j = 0; j < CSEC_BYTES; j++) {
        ASSERT_EQ(tmp[j], (rec.response_outer.get() + i * CSEC_BYTES)[j]);
      }
      memset(tmp, 0, CSEC_BYTES);
    } else {
      for (int j = 0; j < CSEC_BYTES; j++) {
        ASSERT_EQ((snd.base_outer.get() + i * CSEC_BYTES)[j], (rec.response_outer.get() + i * CSEC_BYTES)[j]);
      }
    }
  }
}