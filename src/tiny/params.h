#ifndef TINY_TINY_PARAMS_H_
#define TINY_TINY_PARAMS_H_

#include "OTExtension/util/crypto/crypto.h"
#include "prg/random.h"
#include "tiny-util/network.h"
#include "tiny-util/channel.h"

class Params {
public:
  Params(uint8_t* seed, uint64_t num_pre_gates, uint64_t num_pre_inputs, uint64_t num_pre_outputs, std::string ip_address, uint16_t port, uint8_t net_role, zmq::context_t& context, int num_exces, int exec_id = GLOBAL_PARAMS_CHAN, bool optimize_online = FALSE);
  Params(Params& MainParams, uint8_t* seed, uint64_t num_pre_gates, uint64_t num_pre_inputs, uint64_t num_pre_outputs, int exec_id);

  void ComputeCheckFractions();
  void ComputeGateAndAuthNumbers(uint64_t num_pre_gates, uint64_t num_pre_inputs, uint64_t num_pre_outputs);

  crypto crypt;
  PRNG rnd;
  
  uint64_t delta_pos;
  uint64_t num_pre_gates;
  uint64_t num_pre_inputs;
  uint64_t num_pre_outputs;

  int num_bucket;
  int num_auth;
  int num_inp_auth;
  int num_inp_bucket;

  uint64_t num_eval_gates;
  uint64_t num_eval_auths;

  //Check fractions
  int p_a_base;
  int p_a_nom;
  int p_g_base;
  int p_g_nom;

  int p_a;
  int p_g;
  double s_a;
  double s_g;
  double q_a;
  double q_g;
  
  uint64_t Q;
  uint64_t A;
  uint64_t num_garbled_wires;

  uint64_t left_keys_start;
  uint64_t right_keys_start;
  uint64_t out_keys_start;
  uint64_t auth_start;
  uint64_t ot_chosen_start;
  uint64_t lsb_blind_start;
  uint64_t out_lsb_blind_start;


  //Number indicators
  uint64_t num_OT;
  uint64_t num_commits;

  //Non-Protocol related
  int num_cpus;
  int num_execs;
  int exec_id;
  std::string ip_address;
  uint16_t port;
  uint8_t net_role;
  zmq::context_t& context;
  Channel chan;
};

#endif /* TINY_TINY_PARAMS_H_ */