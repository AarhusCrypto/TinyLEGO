#include "tiny/tiny-evaluator.h"

TinyEvaluator::TinyEvaluator(uint8_t seed[], Params& params) :
  Tiny(seed, params),
  raw_gates_data(5 * params.num_eval_gates * CSEC_BYTES),
  raw_auths_data(3 * params.num_eval_auths * CSEC_BYTES),
  eval_gates_ids(params.num_eval_gates),
  eval_auths_ids(params.num_eval_auths),
  commit_seed_OTs(CODEWORD_BITS),
  commit_seed_choices(CODEWORD_BITS),
  commit_receivers(params.num_max_execs),
  commit_shares(params.num_max_execs),
  global_dot_choices(params.num_pre_inputs),
  global_dot_lsb(params.num_pre_inputs),
  global_input_masks(params.num_pre_inputs, CSEC_BYTES),
  commit_shares_out_lsb_blind(params.num_pre_outputs, CODEWORD_BYTES) {

  //Pointers for convenience to the raw data used for storing the produced eval gates and eval auths. Each exec will write to this array in seperate positions and thus filling it completely. Needs to be computed like this to avoid overflow of evla_data_size
  eval_gates.T_G = raw_gates_data.data();
  eval_gates.T_E = eval_gates.T_G + CSEC_BYTES * params.num_eval_gates;
  eval_gates.S_O = eval_gates.T_E + CSEC_BYTES * params.num_eval_gates;
  eval_gates.S_L = eval_gates.S_O + CSEC_BYTES * params.num_eval_gates;
  eval_gates.S_R = eval_gates.S_L + CSEC_BYTES * params.num_eval_gates;

  eval_auths.H_0 = raw_auths_data.data();
  eval_auths.H_1 = eval_auths.H_0 + CSEC_BYTES * params.num_eval_auths;
  eval_auths.S_A = eval_auths.H_1 + CSEC_BYTES * params.num_eval_auths;
}

TinyEvaluator::~TinyEvaluator() {

  chan.close();

  for (int e = 0; e < exec_channels.size(); ++e) {
    exec_channels[e].close();
  }

  end_point.stop();
}

void TinyEvaluator::Connect(std::string ip_address, uint16_t port) {

  end_point.start(ios, ip_address, port, osuCrypto::EpMode::Client, "ep");

  chan = end_point.addChannel("chan", "chan");

  for (int e = 0; e < commit_receivers.size(); ++e) {
    exec_channels.emplace_back(end_point.addChannel("exec_channel_" + std::to_string(e), "exec_channel_" + std::to_string(e)));
  }
}

void TinyEvaluator::Setup() {

  //BaseOTs
  osuCrypto::u64 num_base_OTs = CSEC + SSEC;
  std::vector<std::array<osuCrypto::block, 2>> base_ots(num_base_OTs);

  osuCrypto::NaorPinkas baseOTs;

  baseOTs.send(base_ots, rnd, chan, 1);

  //Extended the base ots and set them for each dot_receiver
  osuCrypto::KosDotExtReceiver temp_dot_reciever;
  temp_dot_reciever.setBaseOts(base_ots);

  for (int exec_id = 0; exec_id < commit_receivers.size(); ++exec_id) {
    dot_receivers.emplace_back(temp_dot_reciever.split());
  }

  //Extended one last time to setup a kos receiver
  osuCrypto::KosOtExtReceiver kos_receiver;
  std::vector<std::array<osuCrypto::block, 2>> currBaseRecvOts(CSEC);
  for (uint32_t i = 0; i < CSEC; ++i) {
    currBaseRecvOts[i][0] = temp_dot_reciever.mGens[i][0].get<osuCrypto::block>();
    currBaseRecvOts[i][1] = temp_dot_reciever.mGens[i][1].get<osuCrypto::block>();
  }
  kos_receiver.setBaseOts(currBaseRecvOts);

  //Run kos OTX and store the resulting NUM_COMMIT_SEED_OT OTs appropriately
  commit_seed_choices.randomize(rnd);

  kos_receiver.receive(commit_seed_choices, commit_seed_OTs, rnd, chan);

  //Setup tmp commit_receiver
  SplitCommitReceiver tmp_receiver;
  tmp_receiver.SetMsgBitSize(CSEC, gen_matrix_path);

  std::vector<osuCrypto::block> string_msgs(CODEWORD_BITS);
  osuCrypto::BitVector string_choices(CODEWORD_BITS);

  for (int i = 0; i < CODEWORD_BITS; ++i) {
    string_msgs[i] = commit_seed_OTs[i];
    string_choices[i] = commit_seed_choices[i];
  }

  tmp_receiver.SetSeedOTs(string_msgs, string_choices);
  tmp_receiver.GetCloneReceivers(commit_receivers.size(), rnd, commit_receivers, exec_rnds);
}

void TinyEvaluator::Preprocess() {

  //Containers for holding pointers to objects used in each exec. For future use
  std::vector<std::future<void>> cnc_execs_finished(params.num_max_execs);
  std::vector<bool> thread_ver_successes(params.num_max_execs, true);

  //Split the number of preprocessed gates and inputs into num_execs executions
  std::vector<int> inputs_from, inputs_to, outputs_from, outputs_to, gates_from, gates_to, gates_inputs_from, gates_inputs_to;
  PartitionBufferFixedNum(inputs_from, inputs_to, params.num_max_execs, params.num_pre_inputs);
  PartitionBufferFixedNum(gates_inputs_from, gates_inputs_to, params.num_max_execs, params.num_pre_inputs / 2);
  PartitionBufferFixedNum(outputs_from, outputs_to, params.num_max_execs, params.num_pre_outputs);
  PartitionBufferFixedNum(gates_from, gates_to, params.num_max_execs, params.num_pre_gates);

  //Concurrency variables used for ensuring that exec_num 0 has received and updated its global_delta commitment. This is needed as all other executions will use the same commitment to global_delta (in exec_num 0).
  std::mutex delta_received_mutex;
  std::condition_variable delta_received_cond_val;
  bool delta_received = false;
  std::tuple<std::mutex&, std::condition_variable&, bool&> delta_checks = make_tuple(std::ref(delta_received_mutex), std::ref(delta_received_cond_val), std::ref(delta_received));

  //Setup up bucketing permutations
  uint8_t bucket_seed[CSEC];
  rnd.get<uint8_t>(bucket_seed, CSEC_BYTES);
  osuCrypto::PRNG bucket_rnd;
  bucket_rnd.SetSeed(load_block(bucket_seed));
  uint8_t bucket_seeds[2 * CSEC_BYTES];
  bucket_rnd.get<uint8_t>(bucket_seeds, 2 * CSEC_BYTES);

  std::vector<uint32_t> permuted_eval_gates_ids(params.num_eval_gates);
  std::vector<uint32_t> permuted_eval_auths_ids(params.num_eval_auths);

  //Initialize the permutation arrays
  std::iota(std::begin(permuted_eval_gates_ids), std::end(permuted_eval_gates_ids), 0);
  std::iota(std::begin(permuted_eval_auths_ids), std::end(permuted_eval_auths_ids), 0);

  PermuteArray(permuted_eval_gates_ids.data(), params.num_eval_gates, bucket_seeds);
  PermuteArray(permuted_eval_auths_ids.data(), params.num_eval_auths, bucket_seeds + CSEC_BYTES);

  bool delta_flipped = false;
  for (int exec_id = 0; exec_id < params.num_max_execs; ++exec_id) {

    //Assign pr. exec variables that are passed along to the current execution thread
    int inp_from = inputs_from[exec_id];
    int inp_to = inputs_to[exec_id];
    int thread_num_pre_inputs = inp_to - inp_from;
    int thread_num_pre_outputs = outputs_to[exec_id] - outputs_from[exec_id];
    int thread_num_pre_gates = gates_to[exec_id] - gates_from[exec_id];

    //Need to create a new params for each execution with the correct num_pre_gates and num_pre_inputs. The exec_id value decides which channel the execution is communicating on, so must match the constructor execution.
    thread_params_vec.emplace_back(params, thread_num_pre_gates, thread_num_pre_inputs, thread_num_pre_outputs, exec_id);

    //Starts the current execution
    cnc_execs_finished[exec_id] = thread_pool.push([this, exec_id, &delta_checks, &delta_flipped, inp_from, inp_to, &permuted_eval_gates_ids, &permuted_eval_auths_ids, &thread_ver_successes] (int id) {

      uint32_t num_ots;
      if (exec_id == 0) {
        num_ots = thread_params_vec[exec_id].num_pre_inputs + SSEC;
      } else {
        num_ots = thread_params_vec[exec_id].num_pre_inputs;
      }

      uint32_t num_commits = thread_params_vec[exec_id].num_garbled_wires + thread_params_vec[exec_id].num_pre_outputs + num_ots;

      BYTEArrayVector input_masks(num_ots, CSEC_BYTES);

      std::vector<osuCrypto::block> msgs(num_ots);
      osuCrypto::BitVector dot_choices(num_ots);
      dot_choices.randomize(exec_rnds[exec_id]);

      dot_receivers[exec_id]->receive(dot_choices, msgs, exec_rnds[exec_id], exec_channels[exec_id]);

      for (int i = 0; i < num_ots; ++i) {
        _mm_storeu_si128((__m128i*) input_masks[i], msgs[i]);
      }

      commit_shares[exec_id] = BYTEArrayVector(num_commits, CODEWORD_BYTES);

      if (!commit_receivers[exec_id].Commit(commit_shares[exec_id], exec_rnds[exec_id], exec_channels[exec_id])) {
        std::cout << "Abort, key commit failed!" << std::endl;
        throw std::runtime_error("Abort, key commit failed!");
      }

      //Run chosen commit
      BYTEArrayVector input_mask_corrections(num_ots, CSEC_BYTES);
      exec_channels[exec_id].recv(input_mask_corrections.data(), input_mask_corrections.size());


      BYTEArrayVector commit_shares_lsb_blind(SSEC, CODEWORD_BYTES);
      if (!commit_receivers[exec_id].Commit(commit_shares_lsb_blind, exec_rnds[exec_id], exec_channels[exec_id], std::numeric_limits<uint32_t>::max(), ALL_RND_LSB_ZERO)) {
        std::cout << "Abort, blind commit failed!" << std::endl;
        throw std::runtime_error("Abort, blind commit failed!");
      }

      osuCrypto::BitVector decommit_lsb(num_ots);
      exec_channels[exec_id].recv(decommit_lsb.data(), decommit_lsb.sizeBytes());


      BYTEArrayVector commit_shares_ot(num_ots, CODEWORD_BYTES);
      for (int i = 0; i < num_ots; ++i) {
        std::copy(commit_shares[exec_id][thread_params_vec[exec_id].ot_chosen_start + i], commit_shares[exec_id][thread_params_vec[exec_id].ot_chosen_start + i + 1], commit_shares_ot[i]);
      }


      if (!commit_receivers[exec_id].BatchDecommitLSB(commit_shares_ot, decommit_lsb, commit_shares_lsb_blind, exec_rnds[exec_id], exec_channels[exec_id])) {
        std::cout << "Abort, blind lsb decommit failed!" << std::endl;
        throw std::runtime_error("Abort, blind lsb decommit failed!");
      }

      //Put global_delta from OTs in delta_pos of commitment scheme. For security reasons we only do this in exec_num 0, as else a malicious sender might send different delta values in each threaded execution. Therefore only exec_num 0 gets a correction and the rest simply update their delta pointer to point into exec_num 0's delta value.
      std::condition_variable& delta_received_cond_val = std::get<1>(delta_checks);
      bool& delta_received = std::get<2>(delta_checks);

      if (exec_id == 0) {
        uint8_t correction_commit_delta[CODEWORD_BYTES + 1];
        exec_channels[exec_id].recv(correction_commit_delta, CODEWORD_BYTES + 1);

        for (int i = 0; i < CODEWORD_BYTES; ++i) {
          commit_shares[exec_id][thread_params_vec[exec_id].delta_pos][i] ^= (correction_commit_delta[i] & commit_seed_choices.data()[i]);
        }

        delta_received = true;
        delta_flipped = correction_commit_delta[CODEWORD_BYTES];
        delta_received_cond_val.notify_all();

      } else {

        std::mutex& delta_received_mutex = std::get<0>(delta_checks);
        std::unique_lock<std::mutex> lock(delta_received_mutex);
        while (!delta_received) {
          delta_received_cond_val.wait(lock);
        }

        std::copy(commit_shares[0][thread_params_vec[0].delta_pos],
                  commit_shares[0][thread_params_vec[0].delta_pos + 1],
                  commit_shares[exec_id][thread_params_vec[exec_id].delta_pos]);
      }

      if (delta_flipped) {
        for (int i = 0; i < num_ots; ++i) {
          XORBit(127, dot_choices[i], input_masks[i]);
        }
      }

      uint32_t prev_inputs = 0;
      for (int i = 0; i < exec_id; ++i) {
        prev_inputs += thread_params_vec[i].num_pre_inputs;
      }

      for (int i = 0; i < thread_params_vec[exec_id].num_pre_inputs; ++i) {

        XOR_128(global_input_masks[prev_inputs + i], input_mask_corrections[i], input_masks[i]); // turns input_mask_corrections[i] into committed value
        global_dot_choices[prev_inputs + i] = dot_choices[i];
        global_dot_lsb[prev_inputs + i] = decommit_lsb[i];
      }

      //////////////////////////////////CNC////////////////////////////////////
      if (exec_id == 0) {
        //Send own values to sender
        BYTEArrayVector cnc_ot_values(SSEC, CSEC_BYTES);
        osuCrypto::BitVector ot_delta_cnc_choices(SSEC);

        std::copy(input_masks[thread_params_vec[exec_id].num_pre_inputs], input_masks[num_ots], cnc_ot_values.data());

        for (int i = 0; i < SSEC; ++i) {
          ot_delta_cnc_choices[i] = dot_choices[thread_params_vec[exec_id].num_pre_inputs + i];
        }

        exec_channels[exec_id].send(cnc_ot_values.data(), cnc_ot_values.size());
        exec_channels[exec_id].send(ot_delta_cnc_choices.data(), ot_delta_cnc_choices.sizeBytes());

        //Compute decommit shares
        BYTEArrayVector chosen_decommit_shares(SSEC, CODEWORD_BYTES);

        for (int i = 0; i < SSEC; ++i) {
          std::copy(commit_shares[exec_id][thread_params_vec[exec_id].ot_chosen_start + thread_params_vec[exec_id].num_pre_inputs + i], commit_shares[exec_id][thread_params_vec[exec_id].ot_chosen_start + thread_params_vec[exec_id].num_pre_inputs + i + 1], chosen_decommit_shares[i]);

          if (ot_delta_cnc_choices[i]) {
            XOR_CodeWords(chosen_decommit_shares[i], commit_shares[exec_id][thread_params_vec[exec_id].delta_pos]);
          }
        }

        //Receive decommits
        BYTEArrayVector decomitted_values(SSEC, CSEC_BYTES);
        if (!commit_receivers[exec_id].Decommit(chosen_decommit_shares, decomitted_values, exec_channels[exec_id])) {
          thread_ver_successes[exec_id] = false;
          std::cout << "Sender decommit fail in OT CNC!" << std::endl;
          throw std::runtime_error("Sender decommit fail in OT CNC!");
        }

        //Apply the corrections
        uint8_t chosen_decommit_val[CSEC_BYTES];
        for (int i = 0; i < SSEC; ++i) {
          XOR_128(chosen_decommit_val, decomitted_values[i], input_mask_corrections[thread_params_vec[exec_id].num_pre_inputs + i]);

          //Check if they match known value
          if (!std::equal(input_masks[thread_params_vec[exec_id].num_pre_inputs + i], input_masks[thread_params_vec[exec_id].num_pre_inputs + i + 1], chosen_decommit_val)) {
            thread_ver_successes[exec_id] = false;
            std::cout << "Sender cheating in OT CNC. Decomitted to wrong values. Did not commit to Delta!" << std::endl;
            throw std::runtime_error("Sender cheating in OT CNC. Decomitted to wrong values. Did not commit to Delta!");
          }
        }
      }
      //////////////////////////////////CNC////////////////////////////////////


      //================== Preprocess Output Blind Values =====================
      uint32_t prev_outputs = 0;
      for (int i = 0; i < exec_id; ++i) {
        prev_outputs += thread_params_vec[exec_id].num_pre_outputs;
      }

      BYTEArrayVector exec_commit_shares_out_lsb_blind(thread_params_vec[exec_id].num_pre_outputs, CODEWORD_BYTES);

      if (!commit_receivers[exec_id].Commit(exec_commit_shares_out_lsb_blind, exec_rnds[exec_id], exec_channels[exec_id], std::numeric_limits<uint32_t>::max(), ALL_RND_LSB_ZERO)) {
        std::cout << "out_lsb_blind commit failed" << std::endl;
        throw std::runtime_error("out_lsb_blind commit failed");
      }

      for (int i = 0; i < thread_params_vec[exec_id].num_pre_outputs; ++i) {
        std::copy(exec_commit_shares_out_lsb_blind[i],
                  exec_commit_shares_out_lsb_blind[i + 1],
                  commit_shares_out_lsb_blind[prev_outputs + i]);
      }
      //================== Preprocess Output Blind Values =====================

      //==========================Receive Gates===============================
      //Sample the seed used to determine all CNC challenges
      uint8_t cnc_seed[CSEC_BYTES];
      exec_rnds[exec_id].get<uint8_t>(cnc_seed, CSEC_BYTES);

      //Receive all garbling data. When received we send the CNC challenge seed
      std::vector<uint8_t> raw_garbling_data(3 * thread_params_vec[exec_id].Q * CSEC_BYTES + 2 * thread_params_vec[exec_id].A * CSEC_BYTES);
      exec_channels[exec_id].recv(raw_garbling_data.data(), raw_garbling_data.size());

      exec_channels[exec_id].asyncSendCopy(cnc_seed, CSEC_BYTES);

      //Assign pointers to the garbling data. Doing this relatively for clarity
      HalfGates gates_data;
      gates_data.T_G = raw_garbling_data.data();
      gates_data.T_E = gates_data.T_G + thread_params_vec[exec_id].Q * CSEC_BYTES;
      gates_data.S_O = gates_data.T_E + thread_params_vec[exec_id].Q * CSEC_BYTES;

      Auths auths_data;
      auths_data.H_0 = gates_data.S_O + thread_params_vec[exec_id].Q * CSEC_BYTES;
      auths_data.H_1 = auths_data.H_0 + thread_params_vec[exec_id].A * CSEC_BYTES;

      //========================Run Cut-and-Choose=============================

      //Sample check gates and check auths along with the challenge inputs to these. SampleChallenges populates all these variables
      osuCrypto::PRNG cnc_rand;
      cnc_rand.SetSeed(load_block(cnc_seed));

      osuCrypto::BitVector cnc_check_gates(thread_params_vec[exec_id].Q);
      osuCrypto::BitVector cnc_check_auths(thread_params_vec[exec_id].A);
      WeightedRandomString(cnc_check_gates.data(), thread_params_vec[exec_id].p_g, cnc_check_gates.sizeBytes(), cnc_rand);
      WeightedRandomString(cnc_check_auths.data(), thread_params_vec[exec_id].p_a, cnc_check_auths.sizeBytes(), cnc_rand);

      int num_check_gates = countSetBits(cnc_check_gates.data(), 0, thread_params_vec[exec_id].Q - 1);;
      int num_check_auths = countSetBits(cnc_check_auths.data(), 0, thread_params_vec[exec_id].A - 1);;

      osuCrypto::BitVector left_cnc_input(num_check_gates);
      osuCrypto::BitVector right_cnc_input(num_check_gates);
      osuCrypto::BitVector out_cnc_input(num_check_gates);
      osuCrypto::BitVector auth_cnc_input(num_check_auths);

      cnc_rand.get<uint8_t>(left_cnc_input.data(), left_cnc_input.sizeBytes());
      cnc_rand.get<uint8_t>(right_cnc_input.data(), right_cnc_input.sizeBytes());
      for (int i = 0; i < num_check_gates; ++i) {
        out_cnc_input[i] = left_cnc_input[i] & right_cnc_input[i];
      }
      cnc_rand.get<uint8_t>(auth_cnc_input.data(), auth_cnc_input.sizeBytes());
      
      //Construct the CNC check shares to be used for later decommit verification.
      int num_checks = 3 * num_check_gates + num_check_auths;
      BYTEArrayVector cnc_computed_shares(num_checks, CODEWORD_BYTES);

      int current_check_auth_num = 0;
      int current_eval_auth_num = 0;
      std::vector<uint32_t> check_auth_ids(num_check_auths);

      //Development dirty check. Can be deleted once we have computed the correct slack bounds in params.
      bool filled_eval_auths = false;
      for (uint32_t i = 0; i < thread_params_vec[exec_id].A; ++i) {
        if (cnc_check_auths[i]) {
          check_auth_ids[current_check_auth_num] = exec_id * (thread_params_vec[exec_id].Q + thread_params_vec[exec_id].A) + thread_params_vec[exec_id].auth_start + i;
          std::copy(commit_shares[exec_id][thread_params_vec[exec_id].auth_start + i], commit_shares[exec_id][thread_params_vec[exec_id].auth_start + i + 1], cnc_computed_shares[current_check_auth_num]);
          if (auth_cnc_input[current_check_auth_num]) {
            XOR_CodeWords(cnc_computed_shares[current_check_auth_num], commit_shares[exec_id][thread_params_vec[exec_id].delta_pos]);
          }
          ++current_check_auth_num;
        } else {
          //Only write to num_eval_auths if it is not yet filled up. Might be wasteful, but easier to handle
          if (current_eval_auth_num < thread_params_vec[exec_id].num_eval_auths) {

            uint32_t target_pos = permuted_eval_auths_ids[thread_params_vec[exec_id].num_eval_auths * exec_id + current_eval_auth_num];

            std::copy(auths_data.H_0 + i * CSEC_BYTES, auths_data.H_0 + i * CSEC_BYTES + CSEC_BYTES, eval_auths.H_0 + target_pos * CSEC_BYTES);
            std::copy(auths_data.H_1 + i * CSEC_BYTES, auths_data.H_1 + i * CSEC_BYTES + CSEC_BYTES, eval_auths.H_1 + target_pos * CSEC_BYTES);

            //Write the actual auth ID to eval_gates_ids in target_pos, which is determined by permuted_eval_auths_ids
            eval_auths_ids[target_pos] = exec_id * (thread_params_vec[exec_id].Q + thread_params_vec[exec_id].A) + thread_params_vec[exec_id].auth_start + i;
          } else {
            //Development dirty check. Can be deleted once we have computed the correct slack bounds in params.
            filled_eval_auths = true;
          }

          ++current_eval_auth_num;
        }
      }

      //Now for the gates
      std::vector<uint32_t> check_gate_ids(num_check_gates);

      bool filled_eval_gates = false;
      int current_check_gate_num = 0;
      int current_eval_gate_num = 0;
      for (uint32_t i = 0; i < thread_params_vec[exec_id].Q; ++i) {
        if (cnc_check_gates[i]) {
          check_gate_ids[current_check_gate_num] = exec_id * (thread_params_vec[exec_id].Q + thread_params_vec[exec_id].A) + thread_params_vec[exec_id].out_keys_start + i;

          //Left
          std::copy(commit_shares[exec_id][thread_params_vec[exec_id].left_keys_start + i], commit_shares[exec_id][thread_params_vec[exec_id].left_keys_start + i + 1], cnc_computed_shares[num_check_auths + current_check_gate_num]);
          if (left_cnc_input[current_check_gate_num]) {
            XOR_CodeWords(cnc_computed_shares[num_check_auths + current_check_gate_num], commit_shares[exec_id][thread_params_vec[exec_id].delta_pos]);
          }
          //Right
          std::copy(commit_shares[exec_id][thread_params_vec[exec_id].right_keys_start + i], commit_shares[exec_id][thread_params_vec[exec_id].right_keys_start + i + 1], cnc_computed_shares[num_check_auths + num_check_gates + current_check_gate_num]);
          if (right_cnc_input[current_check_gate_num]) {
            XOR_CodeWords(cnc_computed_shares[num_check_auths + num_check_gates + current_check_gate_num], commit_shares[exec_id][thread_params_vec[exec_id].delta_pos]);
          }
          //Out
          std::copy(commit_shares[exec_id][thread_params_vec[exec_id].out_keys_start + i], commit_shares[exec_id][thread_params_vec[exec_id].out_keys_start + i + 1], cnc_computed_shares[num_check_auths + 2 * num_check_gates + current_check_gate_num]);

          if (out_cnc_input[current_check_gate_num]) {
            XOR_CodeWords(cnc_computed_shares[num_check_auths + 2 * num_check_gates + current_check_gate_num], commit_shares[exec_id][thread_params_vec[exec_id].delta_pos]);
          }
          ++current_check_gate_num;
        } else {
          //Only write to num_eval_gates if it is not yet filled up. Might be wasteful, but easier to handle
          if (current_eval_gate_num < thread_params_vec[exec_id].num_eval_gates) {

            int target_pos = permuted_eval_gates_ids[thread_params_vec[exec_id].num_eval_gates * exec_id + current_eval_gate_num];
            std::copy(gates_data.T_G + i * CSEC_BYTES, gates_data.T_G + i * CSEC_BYTES + CSEC_BYTES, eval_gates.T_G + target_pos * CSEC_BYTES);
            std::copy(gates_data.T_E + i * CSEC_BYTES, gates_data.T_E + i * CSEC_BYTES + CSEC_BYTES, eval_gates.T_E + target_pos * CSEC_BYTES);
            std::copy(gates_data.S_O + i * CSEC_BYTES, gates_data.S_O + i * CSEC_BYTES + CSEC_BYTES, eval_gates.S_O + target_pos * CSEC_BYTES);

            //Write the actual gate ID to eval_gates_ids in target_pos, which is determined by permuted_eval_gates_ids
            eval_gates_ids[target_pos] = exec_id * (thread_params_vec[exec_id].Q + thread_params_vec[exec_id].A) + thread_params_vec[exec_id].out_keys_start + i;
          } else {
            //Development dirty check. Can be deleted once we have computed the correct slack bounds in params.
            filled_eval_gates = true;
          }

          ++current_eval_gate_num;
        }
      }

      //Development dirty check. Can be deleted once we have computed the correct slack bounds in params.
      if (!filled_eval_auths) {
        std::cout << "Exec_num: " << exec_id << " did not fill eval_auths" << std::endl;
        thread_ver_successes[exec_id] = false;
      }
      if (!filled_eval_gates) {
        std::cout << "Exec_num: " << exec_id << " did not fill eval_gates" << std::endl;
        thread_ver_successes[exec_id] = false;
      }

      //Receive the 2 * num_check_gates + num_check_auths CNC keys
      int num_check_keys_sent = 2 * num_check_gates + num_check_auths;
      BYTEArrayVector all_cnc_keys(num_checks, CSEC_BYTES);

      exec_channels[exec_id].recv(all_cnc_keys.data(), num_check_keys_sent * CSEC_BYTES);

      GarblingHandler gh(thread_params_vec[exec_id]);
      gh.OutputShiftEvaluateGates(gates_data, 0, all_cnc_keys[num_check_auths], all_cnc_keys[num_check_auths + num_check_gates],
                                  all_cnc_keys[num_check_auths + 2 * num_check_gates],
                                  check_gate_ids.data(), num_check_gates, exec_id * (thread_params_vec[exec_id].Q + thread_params_vec[exec_id].A));

      //Verify the received authenticators
      if (!gh.VerifyAuths(auths_data, 0, all_cnc_keys.data(), check_auth_ids.data(), num_check_auths, exec_id * (thread_params_vec[exec_id].Q + thread_params_vec[exec_id].A))) {
        std::cout << "Auth eval failure!" << std::endl;
        thread_ver_successes[exec_id] = false;
      }

      //Start decommit phase using the above-created indices
      if (!commit_receivers[exec_id].BatchDecommit(cnc_computed_shares, all_cnc_keys, exec_rnds[exec_id], exec_channels[exec_id], true)) {
        std::cout << exec_id << std::endl;
        std::cout << "Wrong keys sent!" << std::endl;
        thread_ver_successes[exec_id] = false;
      }
    });
  }

//Wait for all CNC executions to finish
  for (std::future<void>& r : cnc_execs_finished) {
    r.wait();
  }

  //Send the bucketing info
  chan.send(bucket_seed, CSEC_BYTES);

//Setup maps from eval_gates and eval_auths to commit_block and inner block commit index. Needed to construct decommits that span all executions
  IDMap eval_gates_to_blocks(eval_gates_ids, thread_params_vec[0].Q + thread_params_vec[0].A, thread_params_vec[0].out_keys_start);
  IDMap eval_auths_to_blocks(eval_auths_ids, thread_params_vec[0].Q + thread_params_vec[0].A, thread_params_vec[0].auth_start);

//Starts params.num_max_execs parallel executions for preprocessing solderings. We reuse much of the execution specific information from the last parallel executions
  std::vector<std::future<void>> pre_soldering_execs_finished(params.num_max_execs);
  for (int exec_id = 0; exec_id < params.num_max_execs; ++exec_id) {
    int inp_from = inputs_from[exec_id];
    int inp_to = inputs_to[exec_id];
    int ga_inp_from = gates_inputs_from[exec_id];
    int ga_inp_to = gates_inputs_to[exec_id];
    int ga_from = gates_from[exec_id];
    int ga_to = gates_to[exec_id];

    pre_soldering_execs_finished[exec_id] = thread_pool.push([this, exec_id, inp_from, inp_to, ga_inp_from, ga_inp_to, ga_from, ga_to, &eval_gates_to_blocks, &eval_auths_to_blocks, &thread_ver_successes] (int id) {

      int num_gates = thread_params_vec[exec_id].num_pre_gates;
      int num_inputs = thread_params_vec[exec_id].num_pre_inputs;

      //Set some often used variables
      int num_gate_solderings = num_gates * (params.num_bucket - 1) + (num_inputs / 2) * (params.num_inp_bucket - 1);
      int num_auth_solderings = num_gates * params.num_auth;
      int num_inp_auth_solderings = num_inputs * (params.num_inp_auth - 1);
      int num_pre_solderings = 3 * num_gate_solderings + num_auth_solderings + num_inp_auth_solderings;

      //Receive all preprocessed soldering data and point into this for convenience
      BYTEArrayVector decommited_pre_solderings(num_pre_solderings, CSEC_BYTES);
      exec_channels[exec_id].recv(decommited_pre_solderings.data(), decommited_pre_solderings.size());

      // Apply the actual solderings and store the indices
      int curr_head_gate_pos, curr_head_block, curr_head_idx, curr_gate_pos, curr_gate_block, curr_gate_idx, curr_auth_pos, curr_auth_block, curr_auth_idx, curr_head_inp_auth_pos, curr_head_inp_auth_block, curr_head_inp_auth_idx, curr_inp_auth_pos;
      int solder_gate_pos = 0;
      int solder_auth_pos = 0;
      int solder_inp_auth_pos = 0;
      BYTEArrayVector presolder_computed_shares(num_pre_solderings, CODEWORD_BYTES);

      //We first loop over all head gates
      for (int i = ga_from; i < ga_to; ++i) {
        //Our gate-eval procedure always applies left/right input solderings, so the head gates need to have all-zero left/right solderings. The output soldering is already set, so we ignore this
        curr_head_gate_pos = i * params.num_bucket;
        eval_gates_to_blocks.GetExecIDAndIndex(curr_head_gate_pos, curr_head_block, curr_head_idx);

        //Then all of the gates in this head_gate's bucket (notice j starts at 1).
        for (int j = 1; j < params.num_bucket; ++j) {
          curr_gate_pos = curr_head_gate_pos + j;
          eval_gates_to_blocks.GetExecIDAndIndex(curr_gate_pos, curr_gate_block, curr_gate_idx);

          //Left soldering
          std::copy(decommited_pre_solderings[solder_gate_pos], decommited_pre_solderings[solder_gate_pos + 1], eval_gates.S_L + curr_gate_pos * CSEC_BYTES);

          //Decommit shares
          std::copy(commit_shares[curr_gate_block][thread_params_vec[exec_id].left_keys_start + curr_gate_idx], commit_shares[curr_gate_block][thread_params_vec[exec_id].left_keys_start + curr_gate_idx + 1], presolder_computed_shares[solder_gate_pos]);
          XOR_CodeWords(presolder_computed_shares[solder_gate_pos], commit_shares[curr_head_block][thread_params_vec[exec_id].left_keys_start + curr_head_idx]);

          //Right soldering
          std::copy(decommited_pre_solderings[num_gate_solderings + solder_gate_pos], decommited_pre_solderings[num_gate_solderings + solder_gate_pos + 1], eval_gates.S_R + curr_gate_pos * CSEC_BYTES);

          //Decommit shares
          std::copy(commit_shares[curr_gate_block][thread_params_vec[exec_id].right_keys_start + curr_gate_idx], commit_shares[curr_gate_block][thread_params_vec[exec_id].right_keys_start + curr_gate_idx + 1], presolder_computed_shares[num_gate_solderings + solder_gate_pos]);
          XOR_CodeWords(presolder_computed_shares[num_gate_solderings + solder_gate_pos], commit_shares[curr_head_block][thread_params_vec[exec_id].right_keys_start + curr_head_idx]);

          //Out soldering
          XOR_128(eval_gates.S_O + curr_gate_pos * CSEC_BYTES, decommited_pre_solderings[2 * num_gate_solderings  + solder_gate_pos]);

          //Decommit shares
          std::copy(commit_shares[curr_gate_block][thread_params_vec[exec_id].out_keys_start + curr_gate_idx], commit_shares[curr_gate_block][thread_params_vec[exec_id].out_keys_start + curr_gate_idx + 1], presolder_computed_shares[2 * num_gate_solderings + solder_gate_pos]);
          XOR_CodeWords(presolder_computed_shares[2 * num_gate_solderings + solder_gate_pos], commit_shares[curr_head_block][thread_params_vec[exec_id].out_keys_start + curr_head_idx]);

          ++solder_gate_pos;
        }

        //Then all of the authenticators attached to this head_gate. Here j starts at 0 since all auths are attached to head_gate's output wire
        for (int j = 0; j < params.num_auth; ++j) {
          curr_auth_pos = i * params.num_auth + j;
          eval_auths_to_blocks.GetExecIDAndIndex(curr_auth_pos, curr_auth_block, curr_auth_idx);

          //Inp_auth soldering
          std::copy(decommited_pre_solderings[3 * num_gate_solderings + solder_auth_pos], decommited_pre_solderings[3 * num_gate_solderings + solder_auth_pos + 1], eval_auths.S_A + curr_auth_pos * CSEC_BYTES);

          //Decommit shares
          std::copy(commit_shares[curr_auth_block][thread_params_vec[exec_id].auth_start + curr_auth_idx], commit_shares[curr_auth_block][thread_params_vec[exec_id].auth_start + curr_auth_idx + 1], presolder_computed_shares[3 * num_gate_solderings + solder_auth_pos]);
          XOR_CodeWords(presolder_computed_shares[3 * num_gate_solderings + solder_auth_pos], commit_shares[curr_head_block][thread_params_vec[exec_id].out_keys_start + curr_head_idx]);

          ++solder_auth_pos;
        }
      }

      //Same as above for gates, but for input buckets
      for (int i = ga_inp_from; i < ga_inp_to; ++i) {
        curr_head_gate_pos = params.num_pre_gates * params.num_bucket + i * params.num_inp_bucket;
        eval_gates_to_blocks.GetExecIDAndIndex(curr_head_gate_pos, curr_head_block, curr_head_idx);

        //Then all of the gates in this head_gate's bucket (notice j starts at 1).
        for (int j = 1; j < params.num_inp_bucket; ++j) {
          curr_gate_pos = curr_head_gate_pos + j;
          eval_gates_to_blocks.GetExecIDAndIndex(curr_gate_pos, curr_gate_block, curr_gate_idx);

          //Left soldering
          std::copy(decommited_pre_solderings[solder_gate_pos], decommited_pre_solderings[solder_gate_pos + 1], eval_gates.S_L + curr_gate_pos * CSEC_BYTES);

          //Decommit shares
          std::copy(commit_shares[curr_gate_block][thread_params_vec[exec_id].left_keys_start + curr_gate_idx], commit_shares[curr_gate_block][thread_params_vec[exec_id].left_keys_start + curr_gate_idx + 1], presolder_computed_shares[solder_gate_pos]);
          XOR_CodeWords(presolder_computed_shares[solder_gate_pos], commit_shares[curr_head_block][thread_params_vec[exec_id].left_keys_start + curr_head_idx]);

          //Right soldering
          std::copy(decommited_pre_solderings[num_gate_solderings + solder_gate_pos], decommited_pre_solderings[num_gate_solderings + solder_gate_pos + 1], eval_gates.S_R + curr_gate_pos * CSEC_BYTES);

          //Decommit shares
          std::copy(commit_shares[curr_gate_block][thread_params_vec[exec_id].right_keys_start + curr_gate_idx], commit_shares[curr_gate_block][thread_params_vec[exec_id].right_keys_start + curr_gate_idx + 1], presolder_computed_shares[num_gate_solderings + solder_gate_pos]);
          XOR_CodeWords(presolder_computed_shares[num_gate_solderings + solder_gate_pos], commit_shares[curr_head_block][thread_params_vec[exec_id].right_keys_start + curr_head_idx]);

          //Out soldering
          XOR_128(eval_gates.S_O + curr_gate_pos * CSEC_BYTES, decommited_pre_solderings[2 * num_gate_solderings + solder_gate_pos]);

          //Decommit shares
          std::copy(commit_shares[curr_gate_block][thread_params_vec[exec_id].out_keys_start + curr_gate_idx], commit_shares[curr_gate_block][thread_params_vec[exec_id].out_keys_start + curr_gate_idx + 1], presolder_computed_shares[2 * num_gate_solderings + solder_gate_pos]);
          XOR_CodeWords(presolder_computed_shares[2 * num_gate_solderings + solder_gate_pos], commit_shares[curr_head_block][thread_params_vec[exec_id].out_keys_start + curr_head_idx]);

          ++solder_gate_pos;
        }
      }

      //Finally we create the indices for input authentication. This is constructed exactly the same as the bucket_solderings as the principle is the same using a head input auth and then solderings onto this all the num_inp_auth-1 other authenticators
      for (int i = inp_from; i < inp_to; ++i) {
        curr_head_inp_auth_pos = params.num_pre_gates * params.num_auth + i * params.num_inp_auth;
        eval_auths_to_blocks.GetExecIDAndIndex(curr_head_inp_auth_pos, curr_head_inp_auth_block, curr_head_inp_auth_idx);

        for (int j = 1; j < params.num_inp_auth; ++j) {
          curr_inp_auth_pos = curr_head_inp_auth_pos + j;
          eval_auths_to_blocks.GetExecIDAndIndex(curr_inp_auth_pos, curr_auth_block, curr_auth_idx);

          std::copy(decommited_pre_solderings[3 * num_gate_solderings + num_gates * params.num_auth + solder_inp_auth_pos], decommited_pre_solderings[3 * num_gate_solderings + num_gates * params.num_auth + solder_inp_auth_pos + 1], eval_auths.S_A + curr_inp_auth_pos * CSEC_BYTES);

          //Decommit shares
          std::copy(commit_shares[curr_auth_block][thread_params_vec[exec_id].auth_start + curr_auth_idx], commit_shares[curr_auth_block][thread_params_vec[exec_id].auth_start + curr_auth_idx + 1], presolder_computed_shares[3 * num_gate_solderings + num_auth_solderings + solder_inp_auth_pos]);
          XOR_CodeWords(presolder_computed_shares[3 * num_gate_solderings + num_auth_solderings + solder_inp_auth_pos], commit_shares[curr_head_inp_auth_block][thread_params_vec[exec_id].auth_start + curr_head_inp_auth_idx]);

          ++solder_inp_auth_pos;
        }
      }

      //We end by sending the produced solderings and starting the batch decommit procedure which uses commitments from all executions to build the decommits
      if (!commit_receivers[exec_id].BatchDecommit(presolder_computed_shares, decommited_pre_solderings, exec_rnds[exec_id], exec_channels[exec_id], true)) {
        thread_ver_successes[exec_id] = false;
        std::cout << "Preprocessed soldering decommit failed" << std::endl;
      }
    });
  }

//Wait for all executions to finish
  for (std::future<void>& r : pre_soldering_execs_finished) {
    r.wait();
  }

//Check that all executions in both CNC and preprocessed solderings succeeded
  for (bool b : thread_ver_successes) {
    if (!b) {
      throw std::runtime_error("Abort, initial setup failed. Cheating detected");
    }
  }
}

void TinyEvaluator::Offline(std::vector<Circuit*>& circuits, int top_num_execs) {
  std::vector<std::future<void>> top_soldering_execs_finished(top_num_execs);
  std::vector<bool> thread_ver_successes(top_num_execs, true);

  //Split the number of preprocessed gates and inputs into top_num_execs executions
  std::vector<int> circuits_from, circuits_to;
  PartitionBufferFixedNum(circuits_from, circuits_to, top_num_execs, circuits.size());

  int num_gates_needed = 0;
  int num_inp_gates_needed = 0;
  int num_inps_needed = 0;
  int num_outs_needed = 0;
  for (int i = 0; i < circuits.size(); ++i) {
    gates_offset.emplace_back(num_gates_used + num_gates_needed);
    inp_gates_offset.emplace_back(num_inputs_used / 2 + num_inp_gates_needed);
    inputs_offset.emplace_back(num_inputs_used + num_inps_needed);
    outputs_offset.emplace_back(num_outputs_used + num_outs_needed);
    num_gates_needed += circuits[i]->num_and_gates;
    num_inp_gates_needed += circuits[i]->num_const_inp_wires / 2;
    num_inps_needed += circuits[i]->num_inp_wires;
    num_outs_needed += circuits[i]->num_out_wires;
  }

  if ((params.num_pre_gates - num_gates_used) < num_gates_needed) {
    throw std::runtime_error("Not enough garbled gates");
  } else {
    num_gates_used += num_gates_needed;
  }

  //Due to the way we choose our parameters, if there are enough num_inps, then there are also enough for inp_gates.
  if ((params.num_pre_inputs - num_inputs_used) < num_inps_needed) {
    throw std::runtime_error("Not enough input authenticators");
  } else {
    num_inputs_used += num_inps_needed;
  }

  if ((params.num_pre_outputs - num_outputs_used) < num_outs_needed) {
    throw std::runtime_error("Not enough output wires");
  } else {
    num_outputs_used += num_outs_needed;
  }

  //Setup maps from eval_gates and eval_auths to commit_block and inner block commit index. Needed to construct decommits that span all executions
  IDMap eval_gates_to_blocks(eval_gates_ids, thread_params_vec[0].Q + thread_params_vec[0].A, thread_params_vec[0].out_keys_start);
  IDMap eval_auths_to_blocks(eval_auths_ids, thread_params_vec[0].Q + thread_params_vec[0].A, thread_params_vec[0].auth_start);

  for (int exec_id = 0; exec_id < top_num_execs; ++exec_id) {
    int circ_from = circuits_from[exec_id];
    int circ_to = circuits_to[exec_id];

    top_soldering_execs_finished[exec_id] = thread_pool.push([this, exec_id, circ_from, circ_to, num_outs_needed, &circuits, &eval_gates_to_blocks, &eval_auths_to_blocks, &thread_ver_successes] (int id) {

      //Store in global array
      uint32_t global_out_pos = 0;
      for (int c = 0; c < circ_from; ++c) {
        global_out_pos += circuits[c]->num_out_wires;
      }

      for (int c = circ_from; c < circ_to; ++c) {
        Circuit* circuit = circuits[c];
        int gate_offset = gates_offset[c];
        int inp_gate_offset = inp_gates_offset[c];
        int inp_offset = inputs_offset[c];
        int num_top_solderings = 2 * circuit->num_and_gates + circuit->num_const_inp_wires; //the 2* factor cancels out as we can check two inputs pr. input bucket.

        BYTEArrayVector topsolder_computed_shares(num_top_solderings, CODEWORD_BYTES);
        BYTEArrayVector topsolder_computed_shares_tmp(circuit->num_wires, CODEWORD_BYTES);

        int curr_auth_inp_head_pos, curr_inp_head_block, curr_inp_head_idx, curr_head_pos, curr_head_block, curr_head_idx;
        for (int i = 0; i < circuit->num_inp_wires; ++i) {
          curr_auth_inp_head_pos = params.num_pre_gates * thread_params_vec[exec_id].num_auth + (inp_offset + i) * thread_params_vec[exec_id].num_inp_auth;
          eval_auths_to_blocks.GetExecIDAndIndex(curr_auth_inp_head_pos, curr_inp_head_block, curr_inp_head_idx);

          //Build decommit_info
          std::copy(commit_shares[curr_inp_head_block][thread_params_vec[exec_id].auth_start + curr_inp_head_idx], commit_shares[curr_inp_head_block][thread_params_vec[exec_id].auth_start + curr_inp_head_idx + 1], topsolder_computed_shares_tmp[i]);
        }

        int left_inp_start = circuit->num_and_gates;
        int right_inp_start = 2 * circuit->num_and_gates + circuit->num_const_inp_wires / 2;
        for (int i = 0; i < circuit->num_const_inp_wires / 2; ++i) {
          curr_head_pos = params.num_pre_gates * params.num_bucket + (inp_gate_offset + i) * params.num_inp_bucket;
          eval_gates_to_blocks.GetExecIDAndIndex(curr_head_pos, curr_head_block, curr_head_idx);

          std::copy(commit_shares[curr_head_block][thread_params_vec[exec_id].left_keys_start + curr_head_idx], commit_shares[curr_head_block][thread_params_vec[exec_id].left_keys_start + curr_head_idx + 1], topsolder_computed_shares[left_inp_start + i]);
          XOR_CodeWords(topsolder_computed_shares[left_inp_start + i], topsolder_computed_shares_tmp[i]);
          std::copy(commit_shares[curr_head_block][thread_params_vec[exec_id].right_keys_start + curr_head_idx], commit_shares[curr_head_block][thread_params_vec[exec_id].right_keys_start + curr_head_idx + 1], topsolder_computed_shares[right_inp_start + i]);
          XOR_CodeWords(topsolder_computed_shares[right_inp_start + i], topsolder_computed_shares_tmp[circuit->num_const_inp_wires / 2 + i]);
        }

        int curr_and_gate = 0;
        int left_gate_start = 0;
        int right_gate_start = circuit->num_and_gates + circuit->num_const_inp_wires / 2;
        for (int i = 0; i < circuit->num_gates; ++i) {
          Gate g = circuit->gates[i];
          if (g.type == NOT) {
            //Build decommit_info
            std::copy(topsolder_computed_shares_tmp[g.left_wire], topsolder_computed_shares_tmp[g.left_wire + 1], topsolder_computed_shares_tmp[g.out_wire]);

            XOR_CodeWords(topsolder_computed_shares_tmp[g.out_wire], commit_shares[exec_id][thread_params_vec[exec_id].delta_pos]);

          } else if (g.type == XOR) {
            //Build decommit_info
            std::copy(topsolder_computed_shares_tmp[g.left_wire], topsolder_computed_shares_tmp[g.left_wire + 1], topsolder_computed_shares_tmp[g.out_wire]);

            XOR_CodeWords(topsolder_computed_shares_tmp[g.out_wire], topsolder_computed_shares_tmp[g.right_wire]);

          } else if (g.type == AND) {
            curr_head_pos = (gate_offset + curr_and_gate) * thread_params_vec[exec_id].num_bucket;
            eval_gates_to_blocks.GetExecIDAndIndex(curr_head_pos, curr_head_block, curr_head_idx);
            //Build decommit_info
            std::copy(commit_shares[curr_head_block][thread_params_vec[exec_id].out_keys_start + curr_head_idx], commit_shares[curr_head_block][thread_params_vec[exec_id].out_keys_start + curr_head_idx + 1], topsolder_computed_shares_tmp[g.out_wire]);

            std::copy(commit_shares[curr_head_block][thread_params_vec[exec_id].left_keys_start + curr_head_idx], commit_shares[curr_head_block][thread_params_vec[exec_id].left_keys_start + curr_head_idx + 1], topsolder_computed_shares[left_gate_start + curr_and_gate]);

            std::copy(commit_shares[curr_head_block][thread_params_vec[exec_id].right_keys_start + curr_head_idx], commit_shares[curr_head_block][thread_params_vec[exec_id].right_keys_start + curr_head_idx + 1], topsolder_computed_shares[right_gate_start + curr_and_gate]);

            XOR_CodeWords(topsolder_computed_shares[left_gate_start + curr_and_gate], topsolder_computed_shares_tmp[g.left_wire]);

            XOR_CodeWords(topsolder_computed_shares[right_gate_start + curr_and_gate], topsolder_computed_shares_tmp[g.right_wire]);

            ++curr_and_gate;
          }
        }

        BYTEArrayVector topological_solderings(num_top_solderings, CSEC_BYTES);

        exec_channels[exec_id].recv(topological_solderings.data(), topological_solderings.size());

        if (!commit_receivers[exec_id].BatchDecommit(topsolder_computed_shares, topological_solderings, exec_rnds[exec_id], exec_channels[exec_id], true)) {
          thread_ver_successes[exec_id] = false;
          std::cout << "Topological soldering decommit failed" << std::endl;
        }

        for (int i = 0; i < circuit->num_and_gates; ++i) {
          for (int j = 0; j < thread_params_vec[exec_id].num_bucket; ++j) {
            curr_head_pos = (gate_offset + i) * thread_params_vec[exec_id].num_bucket + j;
            XOR_128(eval_gates.S_L + curr_head_pos * CSEC_BYTES, topological_solderings[left_gate_start + i]);

            XOR_128(eval_gates.S_R + curr_head_pos * CSEC_BYTES, topological_solderings[right_gate_start + i]);
          }
        }

        for (int i = 0; i < circuit->num_const_inp_wires / 2; ++i) {
          for (int j = 0; j < thread_params_vec[exec_id].num_inp_bucket; ++j) {
            curr_head_pos = params.num_pre_gates * params.num_bucket + (inp_gate_offset + i) * params.num_inp_bucket;
            XOR_128(eval_gates.S_L + (curr_head_pos + j) * CSEC_BYTES, topological_solderings[left_inp_start + i]);

            XOR_128(eval_gates.S_R + (curr_head_pos + j) * CSEC_BYTES, topological_solderings[right_inp_start + i]);
          }
        }
      }
    });
  }

  for (std::future<void>& r : top_soldering_execs_finished) {
    r.wait();
  }

  for (bool b : thread_ver_successes) {
    if (!b) {
      throw std::runtime_error("Abort, topological soldering failed. Cheating detected");
    }
  }
}

void TinyEvaluator::Online(std::vector<Circuit*>& circuits, std::vector<osuCrypto::BitVector>& inputs, std::vector<osuCrypto::BitVector>& outputs, int eval_num_execs) {

  std::vector<std::future<void>> online_execs_finished(eval_num_execs);
  std::vector<int> circuits_from, circuits_to;

  PartitionBufferFixedNum(circuits_from, circuits_to, eval_num_execs, circuits.size());

  IDMap eval_auths_to_blocks(eval_auths_ids, thread_params_vec[0].Q + thread_params_vec[0].A, thread_params_vec[0].auth_start);
  IDMap eval_gates_to_blocks(eval_gates_ids, thread_params_vec[0].Q + thread_params_vec[0].A, thread_params_vec[0].out_keys_start);

  for (int exec_id = 0; exec_id < eval_num_execs; ++exec_id) {

    int circ_from = circuits_from[exec_id];
    int circ_to = circuits_to[exec_id];
    Params thread_params = thread_params_vec[exec_id];

    online_execs_finished[exec_id] = thread_pool.push([this, &thread_params, exec_id, circ_from, circ_to, &circuits, &inputs, &outputs, &eval_gates_to_blocks, &eval_auths_to_blocks] (int id) {

      Circuit* circuit;
      Gate g;
      __m128i intrin_outs[thread_params_vec[exec_id].num_bucket];
      __m128i intrin_auths[thread_params_vec[exec_id].num_auth];
      int bucket_score[thread_params_vec[exec_id].num_bucket];
      bool all_equal;

      int gate_offset, inp_gate_offset, inp_offset, out_offset, curr_and_gate;
      int curr_auth_inp_head_pos, curr_inp_head_block, curr_inp_head_idx, curr_head_pos, curr_head_block, curr_head_idx, curr_output_pos, curr_output_block, curr_output_idx;

      GarblingHandler gh(thread_params);
      int curr_input, curr_output, ot_commit_block, commit_id, chosen_val_id;
      for (int c = circ_from; c < circ_to; ++c) {
        circuit = circuits[c];
        gate_offset = gates_offset[c];
        inp_gate_offset = inp_gates_offset[c];
        inp_offset = inputs_offset[c];
        out_offset = outputs_offset[c];

        osuCrypto::BitVector e(circuit->num_eval_inp_wires);

        BYTEArrayVector eval_computed_shares_inp(circuit->num_eval_inp_wires, CODEWORD_BYTES);
        BYTEArrayVector eval_computed_shares_out(circuit->num_out_wires, CODEWORD_BYTES);

        BYTEArrayVector const_inp_keys(circuit->num_const_inp_wires, CSEC_BYTES);
        BYTEArrayVector eval_inp_keys(circuit->num_eval_inp_wires, CSEC_BYTES);

        for (int i = 0; i < circuit->num_eval_inp_wires; ++i) {
          curr_input = (inp_offset + i);
          e[i] = inputs[c][i] ^ global_dot_choices[curr_input];
        }

        if (circuit->num_eval_inp_wires != 0) {
          exec_channels[exec_id].asyncSendCopy(e);
        }

        __m128i* intrin_values = new __m128i[circuit->num_wires]; //using raw pointer due to ~25% increase in overall performance. Since the online phase is so computationally efficient even the slightest performance hit is immediately seen. It does not matter in the others phases as they operation on a very different running time scale.

        for (int i = 0; i < circuit->num_eval_inp_wires; ++i) {
          curr_input = (inp_offset + i);
          ot_commit_block = curr_input / thread_params_vec[exec_id].num_pre_inputs;
          commit_id = thread_params_vec[exec_id].ot_chosen_start + curr_input % thread_params_vec[exec_id].num_pre_inputs;

          std::copy(commit_shares[ot_commit_block][commit_id], commit_shares[ot_commit_block][commit_id + 1], eval_computed_shares_inp[i]);

          //Add the input key
          curr_auth_inp_head_pos = params.num_pre_gates * thread_params_vec[exec_id].num_auth + (circuit->num_const_inp_wires + curr_input) * thread_params_vec[exec_id].num_inp_auth;
          eval_auths_to_blocks.GetExecIDAndIndex(curr_auth_inp_head_pos, curr_inp_head_block, curr_inp_head_idx);
          XOR_CodeWords(eval_computed_shares_inp[i], commit_shares[curr_inp_head_block][thread_params_vec[exec_id].auth_start + curr_inp_head_idx]);

          if (e[i]) {
            XOR_CodeWords(eval_computed_shares_inp[i], commit_shares[ot_commit_block][thread_params_vec[exec_id].delta_pos]);
          }
        }

        if (circuit->num_const_inp_wires != 0) {
          exec_channels[exec_id].recv(const_inp_keys.data(), const_inp_keys.size());
        }

        if (circuit->num_eval_inp_wires != 0) {
          if (!commit_receivers[exec_id].Decommit(eval_computed_shares_inp, eval_inp_keys, exec_channels[exec_id])) {
            std::cout << "Abort: Inp Decommit fail! " << std::endl;
            throw std::runtime_error("Abort: Inp Decommit fail!");
          }
        }

        for (int i = 0; i < circuit->num_eval_inp_wires; ++i) {
          curr_input = (inp_offset + i);
          ot_commit_block = curr_input / thread_params_vec[exec_id].num_pre_inputs;
          chosen_val_id = curr_input % thread_params_vec[exec_id].num_pre_inputs;

          uint8_t lsb_zero_key = GetLSB(eval_inp_keys[i]) ^ global_dot_lsb[curr_input] ^ e[i];

          XOR_128(eval_inp_keys[i], global_input_masks[curr_input]);

          if ((GetLSB(eval_inp_keys[i]) ^ lsb_zero_key) != inputs[c][i]) {
            std::cout << "Abort: Wrong eval value keys sent!" << std::endl;
            throw std::runtime_error("Abort: Wrong eval value keys sent!");
          }

          intrin_values[circuit->num_const_inp_wires + i] = _mm_lddqu_si128((__m128i *) (eval_inp_keys[i]));
        }

        //Ensure that constructor sends valid keys. Implemented different than in paper as we here require that ALL input authenticators accept. This has no influence on security as if we abort here it does not leak anything about the evaluators input. Also, the sender knows if the evaluator is going to abort before sending the bad keys, so it leaks nothing.
        for (int i = 0; i < circuit->num_const_inp_wires; ++i) {
          intrin_values[i] = _mm_lddqu_si128((__m128i *) (const_inp_keys[i]));
        }

        for (int i = 0; i < circuit->num_inp_wires; ++i) {
          curr_auth_inp_head_pos = params.num_pre_gates * thread_params_vec[exec_id].num_auth + (inp_offset + i) * thread_params_vec[exec_id].num_inp_auth;

          for (int j = 0; j < thread_params_vec[exec_id].num_inp_auth; ++j) {
            if (!IntrinVerifyAuths(eval_auths, curr_auth_inp_head_pos + j, intrin_values[i], eval_auths_ids[curr_auth_inp_head_pos + j], gh.key_schedule)) {
              std::cout << "Abort: Inp auth fail! " << std::endl;
              throw std::runtime_error("Abort: Inp auth fail!");
            }
          }
        }

/////////////////////////////// DEBUG Input buckets////////////////////////////
#ifdef DEBUG_SOLDERINGS_INP_BUCKETS
        __m128i out_keys[2];
        for (int i = 0; i < circuit->num_const_inp_wires / 2; ++i) {
          int curr_inp_gate_head = params.num_pre_gates * params.num_bucket + (inp_gate_offset + i) * params.num_inp_bucket;
          IntrinShiftEvaluateGates(eval_gates, curr_inp_gate_head, intrin_values[i], intrin_values[circuit->num_const_inp_wires / 2 + i], out_keys[0], eval_gates_ids[curr_inp_gate_head], gh.key_schedule);

          for (int j = 1; j < params.num_inp_bucket; ++j) {
            IntrinShiftEvaluateGates(eval_gates, curr_inp_gate_head + j, intrin_values[i], intrin_values[circuit->num_const_inp_wires / 2 + i], out_keys[1], eval_gates_ids[curr_inp_gate_head + j], gh.key_schedule);
            if (!compare128(out_keys[0], out_keys[1])) {
              std::cout << "input gate fail pos:" << i << std::endl;
            }
          }
        }
#endif
/////////////////////////////// DEBUG Input buckets////////////////////////////

        curr_and_gate = 0;
        for (int i = 0; i < circuit->num_gates; ++i) {
          g = circuit->gates[i];
          if (g.type == NOT) {
            intrin_values[g.out_wire] = intrin_values[g.left_wire];
          } else if (g.type == XOR) {
            intrin_values[g.out_wire] = _mm_xor_si128(intrin_values[g.left_wire], intrin_values[g.right_wire]);
          } else if (g.type == AND) {
            curr_head_pos = (gate_offset + curr_and_gate) * thread_params_vec[exec_id].num_bucket;

            IntrinShiftEvaluateGates(eval_gates, curr_head_pos, intrin_values[g.left_wire], intrin_values[g.right_wire], intrin_values[g.out_wire], eval_gates_ids[curr_head_pos], gh.key_schedule);

            all_equal = true;
            for (int j = 1; j < thread_params_vec[exec_id].num_bucket; ++j) {
              IntrinShiftEvaluateGates(eval_gates, curr_head_pos + j, intrin_values[g.left_wire], intrin_values[g.right_wire], intrin_outs[j - 1], eval_gates_ids[curr_head_pos + j], gh.key_schedule);
              all_equal &= compare128(intrin_values[g.out_wire], intrin_outs[j - 1]);
            }

            if (!all_equal) {
              std::cout << "all outputs not equal for " << curr_and_gate << std::endl;

              std::fill(bucket_score, bucket_score + thread_params_vec[exec_id].num_bucket * sizeof(uint32_t), 0);
              intrin_outs[thread_params_vec[exec_id].num_bucket - 1] = intrin_values[g.out_wire]; //now all keys are in intrin_outs.
              intrin_auths[0] = intrin_outs[0];
              ++bucket_score[0];
              int candidates = 1;
              for (int j = 1; j < thread_params_vec[exec_id].num_bucket; ++j) {
                int comp = 0;
                for (int k = 0; k < candidates; k++) {
                  comp = !compare128(intrin_outs[j], intrin_outs[k]);
                  if (comp == 0) {
                    ++bucket_score[k];
                    break;
                  }
                }
                if (comp != 0) {
                  intrin_auths[candidates] = intrin_outs[j];
                  ++candidates;
                }
              }

              //Check the candidates
              for (int j = 0; j < thread_params_vec[exec_id].num_auth; j++) {
                curr_auth_inp_head_pos = params.num_pre_gates * thread_params_vec[exec_id].num_auth + (inp_offset + curr_and_gate) * thread_params_vec[exec_id].num_inp_auth;

                for (uint32_t k = 0; k < candidates; k++) {
                  int res = IntrinVerifyAuths(eval_auths, curr_auth_inp_head_pos + k, intrin_auths[k], eval_auths_ids[curr_auth_inp_head_pos + k], gh.key_schedule);
                  if (res == 1) { // The key is good
                    ++bucket_score[k];
                  }
                }
              }
              // Find the winner
              int winner_idx = -1;
              int curr_high_score = -1;
              for (int j = 0; j < candidates; j++) {
                if (bucket_score[j] > curr_high_score) {
                  winner_idx = j;
                  curr_high_score = bucket_score[j];
                }
              }

              intrin_values[g.out_wire] = intrin_auths[winner_idx];
            }
            ++curr_and_gate;
          }
        }

        BYTEArrayVector decommit_shares_out(circuit->num_out_wires, CODEWORD_BYTES);

        //Construct output key decommits
        for (int i = 0; i < circuit->num_out_wires; ++i) {
          curr_output = (out_offset + i);

          std::copy(commit_shares_out_lsb_blind[curr_output], commit_shares_out_lsb_blind[curr_output + 1], decommit_shares_out[i]);

          curr_output_pos = (gate_offset + circuit->num_and_gates - circuit->num_out_wires + i) * thread_params_vec[exec_id].num_bucket;
          eval_gates_to_blocks.GetExecIDAndIndex(curr_output_pos, curr_output_block, curr_output_idx);

          XOR_CodeWords(decommit_shares_out[i], commit_shares[curr_output_block][thread_params_vec[exec_id].out_keys_start + curr_output_idx]);
        }

        BYTEArrayVector decommit_lsb(circuit->num_out_wires, CSEC_BYTES);
        commit_receivers[exec_id].Decommit(decommit_shares_out, decommit_lsb, exec_channels[exec_id]);

        for (int i = 0; i < circuit->num_out_wires; ++i) {
          curr_output = (out_offset + i);
          outputs[c][i] = GetLSB(decommit_lsb[i]) ^ GetLSB(intrin_values[circuit->num_wires - circuit->num_out_wires + i]);
        }

        delete[] intrin_values;
      }
    });
  }

  for (std::future<void>& r : online_execs_finished) {
    r.wait();
  }
}
