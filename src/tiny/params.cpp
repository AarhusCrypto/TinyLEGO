#include "tiny/params.h"

Params::Params(uint64_t num_pre_gates, uint64_t num_pre_inputs, uint64_t num_pre_outputs, int num_execs, int exec_id, bool optimize_online) :
// crypt(CSEC, seed),
num_cpus(std::thread::hardware_concurrency()),
num_execs(num_execs),
exec_id(exec_id)
// context(context),
// ip_address(ip_address),
// port(port),
// net_role(net_role),
// chan(ip_address, port + exec_id + 1, port + exec_id + 1 + MAX_TOTAL_PARAMS, net_role, context) {
  {

  // rnd.SetSeed(seed);

  bool found = false;
  for (int i = 0; i < bucket_param_table_size; ++i) {
    if (num_pre_gates > bucket_param_table[i][0]) {
      num_auth = bucket_param_table[i][1];
      num_bucket = bucket_param_table[i][2];
      p_a = bucket_param_table[i][3];
      p_g = bucket_param_table[i][4];
      found = true;
      break; //Do not try worse parameters
    }
  }
  
  if (!found) {
    throw std::runtime_error("Insufficient number of garbled gates requested.");
  }

  if (optimize_online) {
    bool found = false;
    for (int i = 0; i < bucket_param_table_online_size; ++i) {
      if (num_pre_gates > bucket_param_table_online[i][0]) {
        num_auth = bucket_param_table_online[i][1];
        num_bucket = bucket_param_table_online[i][2];
        p_a = bucket_param_table_online[i][3];
        p_g = bucket_param_table_online[i][4];
        found = true;
        break; //Do not try worse parameters
      }
    }
    if (!found) {
      std::cout << "No online param choices available. Defaulting to normal mode." << std::endl;
    }
  }

  num_inp_auth = (2 * num_auth + 1);
  num_inp_bucket = (2 * num_bucket + 1);

  ComputeGateAndAuthNumbers(num_pre_gates, num_pre_inputs, num_pre_outputs);
}

Params::Params(Params& MainParams, uint64_t num_pre_gates, uint64_t num_pre_inputs, uint64_t num_pre_outputs, int exec_id) :
  // crypt(CSEC, seed),
  num_cpus(std::thread::hardware_concurrency()),
  num_execs(MainParams.num_execs),
  exec_id(exec_id)
  // context(MainParams.context),
  // ip_address(MainParams.ip_address),
  // port(MainParams.port),
  // net_role(MainParams.net_role),
  // chan(ip_address, port + exec_id + 1, port + exec_id + 1 + MAX_TOTAL_PARAMS, net_role, context) {
  {
  // rnd.SetSeed(seed);

  num_auth = MainParams.num_auth;
  num_bucket = MainParams.num_bucket;
  num_inp_auth = MainParams.num_inp_auth;
  num_inp_bucket = MainParams.num_inp_bucket;
  p_a = MainParams.p_a;
  p_g = MainParams.p_g;

  ComputeGateAndAuthNumbers(num_pre_gates, num_pre_inputs, num_pre_outputs);
}

void Params::ComputeGateAndAuthNumbers(uint64_t num_pre_gates, uint64_t num_pre_inputs, uint64_t num_pre_outputs) {

  this->num_pre_gates = num_pre_gates;
  this->num_pre_inputs = num_pre_inputs;
  this->num_pre_outputs = num_pre_outputs;

  ComputeCheckFractions(); // Maybe should only sample check gates once, before going into each execution. Else we need bigger slack!

  Q = ceil((num_pre_gates * num_bucket + num_pre_inputs / 2 * num_inp_bucket) * q_g);
  A = ceil((num_pre_gates * num_auth + num_pre_inputs * num_inp_auth) * q_a);

  num_eval_gates = num_pre_gates * num_bucket + num_pre_inputs / 2 * num_inp_bucket;
  num_eval_auths = num_pre_gates * num_auth + num_pre_inputs * num_inp_auth;

  num_garbled_wires = 3 * Q + A + 1;
  left_keys_start = 0;
  right_keys_start = Q;
  out_keys_start = 2 * Q;
  auth_start = 3 * Q;
  delta_pos = 3 * Q + A;
  lsb_blind_start = 3 * Q + A + 1;
  out_lsb_blind_start = 3 * Q + A + 1 + AES_BITS; //num_pre_outputs commits
  ot_chosen_start = 3 * Q + A + 1 + AES_BITS + num_pre_outputs; //num_pre_inputs + s commits

  // num_OT = num_pre_inputs + SSEC;
  num_commits = num_garbled_wires + AES_BITS + num_pre_outputs + num_pre_inputs; //k commitments will be used to blind VerLeak in offline phase.
}

void Params::ComputeCheckFractions() {

  int p_a_nom = 1;
  int p_g_nom = 1;

  float p_a_tmp = (float)p_a_nom / (float)(1 << p_a);
  float p_g_tmp = (float)p_g_nom / (float)(1 << p_g);

  const double _a_g = 2 * num_pre_gates * num_bucket;
  const double _b = SSEC * log(2);

  double _c = SSEC * log(2) * (p_g_tmp - 1);
  double det = _b * _b - 4 * _a_g * _c;
  if (det < 0.0) {
    printf("ERROR: Bad parameters, p_g slack equation has no real roots");
    exit(EXIT_FAILURE);
  } else {
    // two real roots (possibly equal)
    double r1 = ( -_b + sqrt(det)) / (2 * _a_g);
    double r2 = ( -_b - sqrt(det)) / (2 * _a_g);
//        debug("s_g two roots: %f, %f", r1, r2);
    if (r1 > 0 && r2 > 0) {
      s_g = r1 < r2 ? r1 : r2;
    } else if (r1 <= 0 && r2 > 0) {
      s_g = r2;
    } else if (r1 > 0 && r2 <= 0) {
      s_g = r1;
    } else {
      printf("ERROR: Bad parameters, p_g slack equation has no real roots");
      exit(EXIT_FAILURE);
    }
  }

  const double _a_a = 2 * num_pre_inputs * num_auth;
  _c = SSEC * log(2) * (p_a_tmp - 1);
  det = _b * _b - 4 * _a_a * _c;
  if (det < 0.0) {
    printf("ERROR: Bad parameters, p_a slack equation has no real roots");
    exit(EXIT_FAILURE);
  } else {
    // two real roots (possibly equal)
    double r1 = ( -_b + sqrt(det)) / (2 * _a_a);
    double r2 = ( -_b - sqrt(det)) / (2 * _a_a);
//        debug("s_a two roots: %f, %f", r1, r2);
    if (r1 > 0 && r2 > 0) {
      s_a = r1 < r2 ? r1 : r2;
    } else if (r1 <= 0 && r2 > 0) {
      s_a = r2;
    } else if (r1 > 0 && r2 <= 0) {
      s_a = r1;
    } else {
      printf("ERROR: Bad parameters, p_g slack equation has no real roots");
      exit(EXIT_FAILURE);
    }
  }
  q_a = 1 / (1 - p_a_tmp - s_a);
  q_g = 1 / (1 - p_g_tmp - s_g);
}