#include "tiny/tiny-evaluator.h"

TinyEvaluator::TinyEvaluator(uint8_t seed[], Params& params) :
  Tiny(seed, params),
  // ot_rec(params),
  // rot_seeds(std::make_unique<uint8_t[]>(CODEWORD_BITS * CSEC_BYTES)),
  // rot_choices(std::make_unique<uint8_t[]>(BITS_TO_BYTES(CODEWORD_BITS))),
  verleak_bits(std::make_unique<uint8_t[]>(BITS_TO_BYTES(params.num_pre_outputs + params.num_pre_inputs))),
  raw_eval_data(std::make_unique<uint8_t[]>(5 * CSEC_BYTES * params.num_eval_gates + 3 * CSEC_BYTES * params.num_eval_auths)),
  raw_eval_ids(std::make_unique<uint32_t[]>(params.num_eval_gates + params.num_eval_auths)),
  commit_seed_OTs(CODEWORD_BITS),
  commit_seed_choices(CODEWORD_BITS),
  commit_receivers(params.num_execs),
  commit_shares(params.num_execs) {

  //Pointers for convenience to the raw data used for storing the produced eval gates and eval auths. Each exec will write to this array in seperate positions and thus filling it completely. Needs to be computed like this to avoid overflow of evla_data_size
  eval_gates.T_G = raw_eval_data.get();
  eval_gates.T_E = eval_gates.T_G + CSEC_BYTES * params.num_eval_gates;
  eval_gates.S_O = eval_gates.T_E + CSEC_BYTES * params.num_eval_gates;
  eval_gates.S_L = eval_gates.S_O + CSEC_BYTES * params.num_eval_gates;
  eval_gates.S_R = eval_gates.S_L + CSEC_BYTES * params.num_eval_gates;
  eval_gates_ids = raw_eval_ids.get();

  eval_auths.H_0 = eval_gates.S_R + CSEC_BYTES * params.num_eval_gates;
  eval_auths.H_1 = eval_auths.H_0 + CSEC_BYTES * params.num_eval_auths;
  eval_auths.S_A = eval_auths.H_1 + CSEC_BYTES * params.num_eval_auths;
  eval_auths_ids = eval_gates_ids + params.num_eval_gates;
}

TinyEvaluator::~TinyEvaluator() {

  chan->close();

  for (int e = 0; e < exec_channels.size(); ++e) {
    exec_channels[e]->close();
  }

  end_point.stop();
}

void TinyEvaluator::Connect(std::string ip_address, uint16_t port) {

  end_point.start(ios, ip_address, port, false, "ep");

  chan = &end_point.addChannel("chan", "chan");

  for (int e = 0; e < commit_receivers.size(); ++e) {
    exec_channels.emplace_back(&end_point.addChannel("exec_channel_" + std::to_string(e), "exec_channel_" + std::to_string(e)));
  }
}

void TinyEvaluator::Setup() {
  // ot_rec.net.m_vSocket->ResetSndCnt();
  // ot_rec.net.m_vSocket->ResetRcvCnt();
  //============================Run DOT========================================
  auto baseOT_begin = GET_TIME();
  // ot_rec.InitOTReceiver();

  //BaseOTs
  osuCrypto::u64 num_base_OTs = CSEC + SSEC;
  std::vector<std::array<osuCrypto::block, 2>> base_ots(num_base_OTs);

  osuCrypto::NaorPinkas baseOTs;

  baseOTs.send(base_ots, rnd, *chan, 1);

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

  kos_receiver.receive(commit_seed_choices, commit_seed_OTs, rnd, *chan);

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

  auto baseOT_end = GET_TIME();

#ifdef TINY_PRINT
  PRINT_TIME(baseOT_end, baseOT_begin, "BASEOT");
#endif
}

void TinyEvaluator::Preprocess() {

  std::vector<std::vector<std::chrono::duration<long double, std::milli>>> durations(EVAL_NUM_TIMINGS);

  for (std::vector<std::chrono::duration<long double, std::milli>>& duration : durations) {
    duration.resize(params.num_execs);
  }
  // auto setup_begin = GET_TIME();


  // //Run DOT Extension.

  // ot_rec.Receive();




  // //Construct the CODEWORD_BITS ROTs necessary for commitment scheme. Need to read the bits this way as we are not sure rot_start_pos is 8-bit aligned, if it is not we cannot point at the starting bit.
  // rot_start_pos = params.num_OT - CODEWORD_BITS;
  // for (int i = 0; i < CODEWORD_BITS; ++i) {
  //   if (GetBit(rot_start_pos + i, ot_rec.choices_outer.get())) {
  //     SetBit(i, 1, rot_choices.get());
  //   } else {
  //     SetBit(i, 0, rot_choices.get());
  //   }
  // }

  // //Hash away the global correlation and store in rot_seeds
  // for (int i = 0; i < CODEWORD_BITS; ++i) {
  //   params.crypt.hash(rot_seeds.get() + i * CSEC_BYTES, CSEC_BYTES, ot_rec.response_outer.get() + (rot_start_pos + i) * CSEC_BYTES, CSEC_BYTES);
  // }

  // auto dot_end = GET_TIME();

// #ifdef TINY_PRINT
//   PRINT_TIME(dot_end, dot_begin, "DOT");
// #endif
  //=============================Run Commit====================================
  //Containers for holding pointers to objects used in each exec. For future use
  std::vector<std::future<void>> cnc_execs_finished(params.num_execs);
  std::vector<std::unique_ptr<bool>> thread_ver_successes;

  //Split the number of preprocessed gates and inputs into num_execs executions
  std::vector<int> inputs_from, inputs_to, outputs_from, outputs_to, gates_from, gates_to, gates_inputs_from, gates_inputs_to;
  PartitionBufferFixedNum(inputs_from, inputs_to, params.num_execs, params.num_pre_inputs);
  PartitionBufferFixedNum(gates_inputs_from, gates_inputs_to, params.num_execs, params.num_pre_inputs / 2);
  PartitionBufferFixedNum(outputs_from, outputs_to, params.num_execs, params.num_pre_outputs);
  PartitionBufferFixedNum(gates_from, gates_to, params.num_execs, params.num_pre_gates);

  //Concurrency variables used for ensuring that exec_num 0 has received and updated its global_delta commitment. This is needed as all other executions will use the same commitment to global_delta (in exec_num 0).
  std::mutex cout_mutex;
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

  std::unique_ptr<uint32_t[]> permuted_eval_ids_ptr(new uint32_t[params.num_eval_gates + params.num_eval_auths]);
  uint32_t* permuted_eval_gates_ids = permuted_eval_ids_ptr.get();
  uint32_t* permuted_eval_auths_ids = permuted_eval_gates_ids + params.num_eval_gates;

  //Initialize the permutation arrays
  for (uint32_t i = 0; i < params.num_eval_gates; ++i) {
    permuted_eval_gates_ids[i] = i;
  }
  for (uint32_t i = 0; i < params.num_eval_auths; ++i) {
    permuted_eval_auths_ids[i] = i;
  }
  PermuteArray(permuted_eval_gates_ids, params.num_eval_gates, bucket_seeds);
  PermuteArray(permuted_eval_auths_ids, params.num_eval_auths, bucket_seeds + CSEC_BYTES);

  //store last exec_id as this execution performs the Delta-OT CnC step. This step is needed as we need to signal that the last thread execution ensures that the sender indeed committed to the global_delta used in DOT protocol. We do it in the last execution to avoid dealing with any prefix offset for all OT values, ie. we sacrifices the last SSEC OTs.
  int last_exec_id = params.num_execs - 1;
  //Variable to keep track if execution failed
  for (int exec_id = 0; exec_id < params.num_execs; ++exec_id) {

    //Assign pr. exec variables that are passed along to the current execution thread
    int inp_from = inputs_from[exec_id];
    int inp_to = inputs_to[exec_id];
    int thread_num_pre_inputs = inp_to - inp_from;
    int thread_num_pre_outputs = outputs_to[exec_id] - outputs_from[exec_id];
    int thread_num_pre_gates = gates_to[exec_id] - gates_from[exec_id];

    //Need to create a new params for each execution with the correct num_pre_gates and num_pre_inputs. The exec_id value decides which channel the execution is communicating on, so must match the constructor execution.
    thread_params_vec.emplace_back(std::make_unique<Params>(thread_num_pre_gates, thread_num_pre_inputs, thread_num_pre_outputs, params.num_execs, exec_id));

    // Params* thread_params = thread_params_vec[exec_id].get();

    //We store our local state in the containers as we need to access them for future use
    // commit_recs.emplace_back(std::make_unique<CommitReceiver>(*thread_params, rot_seeds.get(), rot_choices.get()));
    // CommitReceiver* commit_rec = commit_recs[exec_id].get();

    thread_ver_successes.emplace_back(std::make_unique<bool>(true));
    bool* ver_success = thread_ver_successes[exec_id].get();

    //The delta_holder is fixed and passed to all executions as it is in exec 0 that we commit to the global_delta and all other executions use this commitment when needed
    // CommitReceiver* delta_holder = commit_recs[0].get();

    //Starts the current execution
    cnc_execs_finished[exec_id] = thread_pool.push([this, exec_id, &cout_mutex, &delta_checks, ver_success, inp_from, inp_to, last_exec_id, permuted_eval_gates_ids, permuted_eval_auths_ids, &durations] (int id) {

      auto dot_begin = GET_TIME();

      uint32_t num_ots;
      if (exec_id == last_exec_id) {
        num_ots = thread_params_vec[exec_id]->num_pre_inputs + SSEC;
        thread_params_vec[exec_id]->num_commits += SSEC;
      } else {
        num_ots = thread_params_vec[exec_id]->num_pre_inputs;
      }

      BYTEArrayVector input_masks(num_ots, CSEC_BYTES);
      BYTEArrayVector input_masks_choices(BITS_TO_BYTES(num_ots), 1);

      std::vector<osuCrypto::block> msgs(num_ots);
      osuCrypto::BitVector dot_choice(num_ots);
      dot_choice.randomize(exec_rnds[exec_id]);

      dot_receivers[exec_id]->receive(dot_choice, msgs, exec_rnds[exec_id], *exec_channels[exec_id]);

      for (int i = 0; i < num_ots; ++i) {
        _mm_storeu_si128((__m128i*) input_masks[i], msgs[i]);
        SetBit(i, dot_choice[i], input_masks_choices.data());
      }

      auto dot_end = GET_TIME();


      auto commit_begin = GET_TIME();

      //If it's the last execution then we commit to s extra OTs as these are to be used for CNC.
      // int num_OT_commits = inp_to - inp_from;
      // if (exec_id == last_exec_id) {
      //   num_OT_commits += SSEC;
      //   thread_params_vec[exec_id]->num_commits += SSEC;
      // }

      // Do chosen commit to all DOT commitments.
      // commit_rec->ChosenCommit(num_OT_commits);

      commit_shares[exec_id] = BYTEArrayVector(thread_params_vec[exec_id]->num_commits, CODEWORD_BYTES);

      if (!commit_receivers[exec_id].Commit(commit_shares[exec_id], exec_rnds[exec_id], *exec_channels[exec_id])) {
        std::cout << "Abort, key commit failed!" << std::endl;
        throw std::runtime_error("Abort, key commit failed!");
      }

      //Run chosen commit
      BYTEArrayVector input_mask_corrections(num_ots, CSEC_BYTES);
      exec_channels[exec_id]->recv(input_mask_corrections.data(), input_mask_corrections.size());

      auto commit_end = GET_TIME();
      durations[EVAL_COMMIT_TIME][exec_id] = commit_end - commit_begin;

      //Put global_delta from OTs in delta_pos of commitment scheme. For security reasons we only do this in exec_num 0, as else a malicious sender might send different delta values in each threaded execution. Therefore only exec_num 0 gets a correction and the rest simply update their delta pointer to point into exec_num 0's delta value.
      std::condition_variable& delta_received_cond_val = std::get<1>(delta_checks);
      bool& delta_received = std::get<2>(delta_checks);
      if (exec_id == 0) {
        uint8_t correction_commit_delta[CODEWORD_BYTES];
        exec_channels[exec_id]->recv(correction_commit_delta, CODEWORD_BYTES);

        for (int i = 0; i < CODEWORD_BYTES; ++i) {
          commit_shares[exec_id][thread_params_vec[exec_id]->delta_pos][i] ^= (correction_commit_delta[i] & commit_seed_choices.data()[i]);
        }

        delta_received = true;
        delta_received_cond_val.notify_all();

      } else {

        std::mutex& delta_received_mutex = std::get<0>(delta_checks);
        std::unique_lock<std::mutex> lock(delta_received_mutex);
        while (!delta_received) {
          delta_received_cond_val.wait(lock);
        }
        std::copy(commit_shares[0][thread_params_vec[0]->delta_pos],
                  commit_shares[0][thread_params_vec[0]->delta_pos + 1],
                  commit_shares[exec_id][thread_params_vec[exec_id]->delta_pos]);
      }

      //////////////////////////////////////CNC////////////////////////////////////
      if (exec_id == last_exec_id) {
        //Send own values to sender
        std::vector<uint8_t> cnc_ot_values(SSEC * CSEC_BYTES + SSEC_BYTES);
        uint8_t* ot_delta_cnc_choices = cnc_ot_values.data() + SSEC * CSEC_BYTES;

        std::copy(input_masks[thread_params_vec[exec_id]->num_pre_inputs], input_masks[num_ots], cnc_ot_values.data());

        for (int i = 0; i < SSEC; ++i) {
          if (GetBit((thread_params_vec[exec_id]->num_pre_inputs + i), input_masks_choices.data())) {
            SetBit(i, 1, ot_delta_cnc_choices);
          } else {
            SetBit(i, 0, ot_delta_cnc_choices);
          }
        }

        exec_channels[exec_id]->asyncSendCopy(cnc_ot_values.data(),  SSEC * CSEC_BYTES + SSEC_BYTES);

        //Compute decommit shares
        BYTEArrayVector chosen_decommit_shares(SSEC, CODEWORD_BYTES);
        
        for (int i = 0; i < SSEC; ++i) {
          std::copy(commit_shares[exec_id][thread_params_vec[exec_id]->ot_chosen_start + thread_params_vec[exec_id]->num_pre_inputs + i], commit_shares[exec_id][thread_params_vec[exec_id]->ot_chosen_start + thread_params_vec[exec_id]->num_pre_inputs + i + 1], chosen_decommit_shares[i]);

          if (GetBit(i, ot_delta_cnc_choices)) {
            XOR_CodeWords(chosen_decommit_shares[i], commit_shares[exec_id][thread_params_vec[exec_id]->delta_pos]);
          }
        }

        //Receive decommits
        BYTEArrayVector decomitted_values(SSEC, CSEC_BYTES);
        if (!commit_receivers[exec_id].Decommit(chosen_decommit_shares, decomitted_values, *exec_channels[exec_id])) {
          *ver_success = false;
          std::cout << "Sender decommit fail in OT CNC!" << std::endl;
          throw std::runtime_error("Sender decommit fail in OT CNC!");
        }

        //Apply the corrections
        uint8_t chosen_decommit_val[CSEC_BYTES];
        for (int i = 0; i < SSEC; ++i) {
          XOR_128(chosen_decommit_val, decomitted_values[i], input_mask_corrections[thread_params_vec[exec_id]->num_pre_inputs + i]);

          //Check if they match known value
          if (!std::equal(input_masks[thread_params_vec[exec_id]->num_pre_inputs + i], input_masks[thread_params_vec[exec_id]->num_pre_inputs + i + 1], chosen_decommit_val)) {
            *ver_success = false;
            std::cout << "Sender cheating in OT CNC. Decomitted to wrong values. Did not commit to Delta!" << std::endl;
            throw std::runtime_error("Sender cheating in OT CNC. Decomitted to wrong values. Did not commit to Delta!");
          }
        }
      }
//////////////////////////////////////CNC////////////////////////////////////

      for (int i = 0; i < thread_params_vec[exec_id]->num_pre_inputs; ++i) {

        XOR_128(input_mask_corrections[i], input_masks[i]); // turns input_mask_corrections[i] into committed value
      }

      // //If in the last execution we do CNC on the global_delta to ensure that this is indeed the global delta used in DOT as well
      // if (exec_id == last_exec_id) {

      //   //Send own OT value for the last s OTs + own OT choice-bit
      //   std::unique_ptr<uint8_t[]> cnc_ot_values(std::make_unique<uint8_t[]>(SSEC * CSEC_BYTES + SSEC_BYTES));
      //   std::copy(ot_rec.response_outer.get() + inp_to * CSEC_BYTES, ot_rec.response_outer.get() + inp_to * CSEC_BYTES + SSEC * CSEC_BYTES, cnc_ot_values.get());
      //   uint8_t* delta_ot_cnc_choices = cnc_ot_values.get() + SSEC * CSEC_BYTES;
      //   for (int i = 0; i < SSEC; ++i) {
      //     if (GetBit((inp_to + i), ot_rec.choices_outer.get())) {
      //       SetBit(i, 1, delta_ot_cnc_choices);
      //     } else {
      //       SetBit(i, 0, delta_ot_cnc_choices);
      //     }
      //   }
      //   thread_params->chan.Send(cnc_ot_values.get(),  SSEC * CSEC_BYTES + SSEC_BYTES);

      //   //Build decommit info
      //   std::unique_ptr<uint8_t[]> OT_decomitted_values(std::make_unique<uint8_t[]>(SSEC * CSEC_BYTES));
      //   std::unique_ptr<uint8_t[]> chosen_computed_shares(std::make_unique<uint8_t[]>(SSEC * CODEWORD_BYTES));

      //   for (int i = 0; i < SSEC; ++i) {
      //     int commit_id = thread_params->ot_chosen_start + num_OT_commits - SSEC + i;
      //     std::copy(commit_rec->commit_shares[commit_id], commit_rec->commit_shares[commit_id] + CODEWORD_BYTES, chosen_computed_shares.get() + i * CODEWORD_BYTES);
      //     if (GetBit(i, delta_ot_cnc_choices)) {
      //       XOR_CodeWords(chosen_computed_shares.get() + i * CODEWORD_BYTES, commit_rec->commit_shares[thread_params->delta_pos]);
      //     }
      //   }

      //   std::vector<uint64_t> chosen_decommit_pos(SSEC);
      //   std::iota(std::begin(chosen_decommit_pos), std::end(chosen_decommit_pos), num_OT_commits - SSEC);

      //   if (commit_rec->ChosenDecommit(chosen_computed_shares.get(), OT_decomitted_values.get(), chosen_decommit_pos, SSEC)) {
      //     for (int i = 0; i < SSEC; ++i) {
      //       uint8_t* own_ot = ot_rec.response_outer.get() + inp_to * CSEC_BYTES + i * CSEC_BYTES;
      //       uint8_t* decom_ot = OT_decomitted_values.get() + i * CSEC_BYTES;
      //       if (!std::equal(own_ot, own_ot + CSEC_BYTES, decom_ot)) {
      //         std::cout << "Sender cheating in OT CNC. Decomitted to wrong values!" << std::endl;
      //         *ver_success = false;
      //       }
      //     }
      //   } else {
      //     std::cout << "Sender cheating in OT CNC. Decommit failed!" << std::endl;
      //     *ver_success = false;
      //   }
      // }

//       //===========================VER_LEAK====================================
//       auto verleak_begin = GET_TIME();
//       int num_verleaks = thread_params->num_pre_outputs + thread_params->num_pre_inputs;

//       std::unique_ptr<uint8_t[]> verleak_computed_shares(std::make_unique<uint8_t[]>((num_verleaks + AES_BITS) * CODEWORD_BYTES + BITS_TO_BYTES(num_verleaks + AES_BITS)));
//       for (int i = 0; i < num_verleaks; ++i) {
//         std::copy(commit_rec->commit_shares[thread_params->out_lsb_blind_start + i], commit_rec->commit_shares[thread_params->out_lsb_blind_start + i] + CODEWORD_BYTES, verleak_computed_shares.get() + i * CODEWORD_BYTES);
//       }

//       for (int i = 0; i < AES_BITS; ++i) {
//         std::copy(commit_rec->commit_shares[thread_params->lsb_blind_start + i], commit_rec->commit_shares[thread_params->lsb_blind_start + i] + CODEWORD_BYTES, verleak_computed_shares.get() + (num_verleaks + i) * CODEWORD_BYTES);
//       }

//       uint8_t* local_verleak_bits = verleak_computed_shares.get() + (num_verleaks + AES_BITS) * CODEWORD_BYTES;
//       thread_params->chan.ReceiveBlocking(local_verleak_bits, BITS_TO_BYTES(num_verleaks + AES_BITS));

//       if (!BatchDecommitLSB(commit_rec, verleak_computed_shares.get(), num_verleaks + AES_BITS, local_verleak_bits)) {
//         std::cout << "VerLeak Failed!" << std::endl;
//         *ver_success = false;
//       }

//       //Update and check the ver_leak bits for the input masking values
//       for (int i = 0; i < thread_params->num_pre_inputs; ++i) {
//         XORBit(thread_params->num_pre_outputs + i, GetLSB(commit_rec->chosen_commit_values.get() + i * CSEC_BYTES), local_verleak_bits);

//         //When done use the last num_pre_inputs bits in a verification step involving outer_response and outer_choice from the OT. This will catch if sender did not commit to the 0-value.
//         if ((GetLSB(ot_rec.response_outer.get() + (inp_from + i) * CSEC_BYTES) ^ GetBit(inp_from + i, ot_rec.choices_outer.get())) != (GetBit(thread_params->num_pre_outputs + i, local_verleak_bits))) {
//           std::cout << "Wrong lsb values sent!" << std::endl;
//           *ver_success = false;
//         }
//       }

//       //Write to global ver_leak array
//       for (int i = 0; i < thread_params->num_pre_outputs; ++i) {
//         SetBit(exec_id * thread_params->num_pre_outputs + i, GetBit(i, local_verleak_bits), verleak_bits.get());
//       }

//       for (int i = 0; i < thread_params->num_pre_inputs; ++i) {
//         SetBit(params.num_pre_outputs + exec_id * thread_params->num_pre_inputs + i, GetBit(thread_params->num_pre_outputs + i, local_verleak_bits), verleak_bits.get());
//       }
//       auto verleak_end = GET_TIME();
//       durations[EVAL_VERLEAK_TIME][exec_id] = verleak_end - verleak_begin;
//       //==========================Receive Gates===============================
//       auto receive_gates_auths_begin = GET_TIME();

//       //Sample the seed used to determine all CNC challenges
//       uint8_t cnc_seed[CSEC_BYTES];
//       thread_params->rnd.GenRnd(cnc_seed, CSEC_BYTES);

//       //Receive all garbling data. When received we send the CNC challenge seed
//       std::unique_ptr<uint8_t[]> raw_garbling_data(std::make_unique<uint8_t[]>(3 * thread_params->Q * CSEC_BYTES + 2 * thread_params->A * CSEC_BYTES));
//       thread_params->chan.ReceiveBlocking(raw_garbling_data.get(), 3 * thread_params->Q * CSEC_BYTES + 2 * thread_params->A * CSEC_BYTES);

//       thread_params->chan.Send(cnc_seed, CSEC_BYTES);

//       //Assign pointers to the garbling data. Doing this relatively for clarity
//       HalfGates gates_data;
//       gates_data.T_G = raw_garbling_data.get();
//       gates_data.T_E = gates_data.T_G + thread_params->Q * CSEC_BYTES;
//       gates_data.S_O = gates_data.T_E + thread_params->Q * CSEC_BYTES;

//       Auths auths_data;
//       auths_data.H_0 = gates_data.S_O + thread_params->Q * CSEC_BYTES;
//       auths_data.H_1 = auths_data.H_0 + thread_params->A * CSEC_BYTES;

//       auto receive_gates_auths_end = GET_TIME();
//       durations[EVAL_RECEIVE_GATES_AUTHS_TIME][exec_id] = receive_gates_auths_end - receive_gates_auths_begin;

//       //========================Run Cut-and-Choose=============================
//       auto cnc_begin = GET_TIME();

//       //Sample check gates and check auths along with the challenge inputs to these. SampleChallenges populates all these variables
//       int num_bytes_gates = BITS_TO_BYTES(thread_params->Q);
//       int num_bytes_auths = BITS_TO_BYTES(thread_params->A);
//       PRNG cnc_rand;
//       cnc_rand.SetSeed(cnc_seed);

//       std::unique_ptr<uint8_t[]> cnc_check_gates(std::make_unique<uint8_t[]>(num_bytes_gates + num_bytes_auths));
//       uint8_t* cnc_check_auths = cnc_check_gates.get() + num_bytes_gates;
//       WeightedRandomString(cnc_check_gates.get(), thread_params->p_g, num_bytes_gates, cnc_rand);
//       WeightedRandomString(cnc_check_auths, thread_params->p_a, num_bytes_auths, cnc_rand);

//       int num_check_gates = countSetBits(cnc_check_gates.get(), 0, thread_params->Q - 1);
//       int num_check_auths = countSetBits(cnc_check_auths, 0, thread_params->A - 1);

//       std::unique_ptr<uint8_t[]> left_cnc_input(std::make_unique<uint8_t[]>(3 * BITS_TO_BYTES(num_check_gates) + BITS_TO_BYTES(num_check_auths)));
//       uint8_t* right_cnc_input = left_cnc_input.get() + BITS_TO_BYTES(num_check_gates);
//       uint8_t* out_cnc_input = right_cnc_input + BITS_TO_BYTES(num_check_gates);
//       uint8_t* auth_cnc_input = out_cnc_input + BITS_TO_BYTES(num_check_gates);

//       cnc_rand.GenRnd(left_cnc_input.get(), BITS_TO_BYTES(num_check_gates));
//       cnc_rand.GenRnd(right_cnc_input, BITS_TO_BYTES(num_check_gates));
//       for (int i = 0; i < BITS_TO_BYTES(num_check_gates); ++i) {
//         out_cnc_input[i] = left_cnc_input[i] & right_cnc_input[i];
//       }
//       cnc_rand.GenRnd(auth_cnc_input, BITS_TO_BYTES(num_check_auths));


//       //Construct the CNC check shares to be used for later decommit verification.
//       int num_checks = 3 * num_check_gates + num_check_auths;
//       std::unique_ptr<uint8_t[]> cnc_computed_shares(std::make_unique<uint8_t[]>(num_checks * CODEWORD_BYTES));

//       int current_check_auth_num = 0;
//       int current_eval_auth_num = 0;
//       std::unique_ptr<uint32_t[]> check_auth_ids(std::make_unique<uint32_t[]>(num_check_auths));

//       //Development dirty check. Can be deleted once we have computed the correct slack bounds in params.
//       bool filled_eval_auths = false;
//       for (uint32_t i = 0; i < thread_params->A; ++i) {
//         if (GetBit(i, cnc_check_auths)) {
//           check_auth_ids[current_check_auth_num] = exec_id * (thread_params->Q + thread_params->A) + thread_params->auth_start + i;
//           std::copy(commit_rec->commit_shares[thread_params->auth_start + i], commit_rec->commit_shares[thread_params->auth_start + i] + CODEWORD_BYTES, cnc_computed_shares.get() + current_check_auth_num * CODEWORD_BYTES);
//           if (GetBit(current_check_auth_num, auth_cnc_input)) {
//             XOR_CodeWords(cnc_computed_shares.get() + current_check_auth_num * CODEWORD_BYTES, commit_rec->commit_shares[thread_params->delta_pos]);
//           }
//           ++current_check_auth_num;
//         } else {
//           //Only write to num_eval_auths if it is not yet filled up. Might be wasteful, but easier to handle
//           if (current_eval_auth_num < thread_params->num_eval_auths) {

//             uint32_t target_pos = permuted_eval_auths_ids[thread_params->num_eval_auths * exec_id + current_eval_auth_num];
//             std::copy(auths_data.H_0 + i * CSEC_BYTES, auths_data.H_0 + i * CSEC_BYTES + CSEC_BYTES, eval_auths.H_0 + target_pos * CSEC_BYTES);
//             std::copy(auths_data.H_1 + i * CSEC_BYTES, auths_data.H_1 + i * CSEC_BYTES + CSEC_BYTES, eval_auths.H_1 + target_pos * CSEC_BYTES);

//             //Write the actual auth ID to eval_gates_ids in target_pos, which is determined by permuted_eval_auths_ids
//             eval_auths_ids[target_pos] = exec_id * (thread_params->Q + thread_params->A) + thread_params->auth_start + i;
//           } else {
//             //Development dirty check. Can be deleted once we have computed the correct slack bounds in params.
//             filled_eval_auths = true;
//           }

//           ++current_eval_auth_num;
//         }
//       }

//       //Now for the gates
//       std::unique_ptr<uint32_t[]> check_gate_ids(std::make_unique<uint32_t[]>(num_check_gates));

//       bool filled_eval_gates = false;
//       int current_check_gate_num = 0;
//       int current_eval_gate_num = 0;
//       for (uint32_t i = 0; i < thread_params->Q; ++i) {
//         if (GetBit(i, cnc_check_gates.get())) {
//           check_gate_ids[current_check_gate_num] = exec_id * (thread_params->Q + thread_params->A) + thread_params->out_keys_start + i;

//           //Left
//           std::copy(commit_rec->commit_shares[thread_params->left_keys_start + i], commit_rec->commit_shares[thread_params->left_keys_start + i] + CODEWORD_BYTES, cnc_computed_shares.get() + (num_check_auths + current_check_gate_num) * CODEWORD_BYTES);
//           if (GetBit(current_check_gate_num, left_cnc_input.get())) {
//             XOR_CodeWords(cnc_computed_shares.get() + (num_check_auths + current_check_gate_num) * CODEWORD_BYTES, commit_rec->commit_shares[thread_params->delta_pos]);
//           }
//           //Right
//           std::copy(commit_rec->commit_shares[thread_params->right_keys_start + i], commit_rec->commit_shares[thread_params->right_keys_start + i] + CODEWORD_BYTES, cnc_computed_shares.get() + (num_check_auths + num_check_gates + current_check_gate_num) * CODEWORD_BYTES);
//           if (GetBit(current_check_gate_num, right_cnc_input)) {
//             XOR_CodeWords(cnc_computed_shares.get() + (num_check_auths + num_check_gates + current_check_gate_num) * CODEWORD_BYTES, commit_rec->commit_shares[thread_params->delta_pos]);
//           }
//           //Out
//           std::copy(commit_rec->commit_shares[thread_params->out_keys_start + i], commit_rec->commit_shares[thread_params->out_keys_start + i] + CODEWORD_BYTES, cnc_computed_shares.get() + (num_check_auths + 2 * num_check_gates + current_check_gate_num) * CODEWORD_BYTES);

//           if (GetBit(current_check_gate_num, out_cnc_input)) {
//             XOR_CodeWords(cnc_computed_shares.get() + (num_check_auths + 2 * num_check_gates + current_check_gate_num) * CODEWORD_BYTES, commit_rec->commit_shares[thread_params->delta_pos]);
//           }
//           ++current_check_gate_num;
//         } else {
//           //Only write to num_eval_gates if it is not yet filled up. Might be wasteful, but easier to handle
//           if (current_eval_gate_num < thread_params->num_eval_gates) {

//             int target_pos = permuted_eval_gates_ids[thread_params->num_eval_gates * exec_id + current_eval_gate_num];
//             std::copy(gates_data.T_G + i * CSEC_BYTES, gates_data.T_G + i * CSEC_BYTES + CSEC_BYTES, eval_gates.T_G + target_pos * CSEC_BYTES);
//             std::copy(gates_data.T_E + i * CSEC_BYTES, gates_data.T_E + i * CSEC_BYTES + CSEC_BYTES, eval_gates.T_E + target_pos * CSEC_BYTES);
//             std::copy(gates_data.S_O + i * CSEC_BYTES, gates_data.S_O + i * CSEC_BYTES + CSEC_BYTES, eval_gates.S_O + target_pos * CSEC_BYTES);

//             //Write the actual gate ID to eval_gates_ids in target_pos, which is determined by permuted_eval_gates_ids
//             eval_gates_ids[target_pos] = exec_id * (thread_params->Q + thread_params->A) + thread_params->out_keys_start + i;
//           } else {
//             //Development dirty check. Can be deleted once we have computed the correct slack bounds in params.
//             filled_eval_gates = true;
//           }

//           ++current_eval_gate_num;
//         }
//       }

//       //Development dirty check. Can be deleted once we have computed the correct slack bounds in params.
//       if (!filled_eval_auths) {
//         std::cout << "Exec_num: " << exec_id << " did not fill eval_auths" << std::endl;
//         *ver_success = false;
//       }
//       if (!filled_eval_gates) {
//         std::cout << "Exec_num: " << exec_id << " did not fill eval_gates" << std::endl;
//         *ver_success = false;
//       }

//       //Receive the 2 * num_check_gates + num_check_auths CNC keys
//       int num_checks_sent = 2 * num_check_gates + num_check_auths;
//       std::unique_ptr<uint8_t[]> all_cnc_keys(std::make_unique<uint8_t[]>(num_checks * CSEC_BYTES));

//       thread_params->chan.ReceiveBlocking(all_cnc_keys.get(), num_checks_sent * CSEC_BYTES);

//       GarblingHandler gh(*thread_params);
//       gh.OutputShiftEvaluateGates(gates_data, 0, all_cnc_keys.get() + num_check_auths * CSEC_BYTES, all_cnc_keys.get() + (num_check_auths + num_check_gates) * CSEC_BYTES, all_cnc_keys.get() + (num_check_auths + 2 * num_check_gates) * CSEC_BYTES, check_gate_ids.get(), num_check_gates, exec_id * (thread_params->Q + thread_params->A));

//       //Verify the received authenticators
//       if (!gh.VerifyAuths(auths_data, 0, all_cnc_keys.get(), check_auth_ids.get(), num_check_auths, exec_id * (thread_params->Q + thread_params->A))) {
//         std::cout << "Auth eval failure!" << std::endl;
//         *ver_success = false;
//       };

//       //Start decommit phase using the above-created indices
//       auto cnc_batch_decommit_begin = GET_TIME();
//       if (!commit_rec->BatchDecommit(cnc_computed_shares.get(), num_checks, all_cnc_keys.get())) { //
//         std::cout << "Wrong keys sent!" << std::endl;
//         *ver_success = false;
//       }
//       auto cnc_end = GET_TIME();
//       durations[EVAL_CNC_TIME][exec_id] = cnc_end - cnc_begin;
    });
  }

//Wait for all CNC executions to finish
  for (std::future<void>& r : cnc_execs_finished) {
    r.wait();
  }

//   //Send the bucketing info
//   params.chan.Send(bucket_seed, CSEC_BYTES);

// //Setup maps from eval_gates and eval_auths to commit_block and inner block commit index. Needed to construct decommits that span all executions
//   IDMap eval_gates_to_blocks(eval_gates_ids, thread_params_vec[0]->Q + thread_params_vec[0]->A, thread_params_vec[0]->out_keys_start);
//   IDMap eval_auths_to_blocks(eval_auths_ids, thread_params_vec[0]->Q + thread_params_vec[0]->A, thread_params_vec[0]->auth_start);

// //Starts params.num_execs parallel executions for preprocessing solderings. We reuse much of the execution specific information from the last parallel executions
//   auto presolder_begin = GET_TIME();
//   std::vector<std::future<void>> pre_soldering_execs_finished(params.num_execs);
//   for (int exec_id = 0; exec_id < params.num_execs; ++exec_id) {
//     int inp_from = inputs_from[exec_id];
//     int inp_to = inputs_to[exec_id];
//     int ga_inp_from = gates_inputs_from[exec_id];
//     int ga_inp_to = gates_inputs_to[exec_id];
//     int ga_from = gates_from[exec_id];
//     int ga_to = gates_to[exec_id];
//     bool* ver_success = thread_ver_successes[exec_id].get();
//     Params* thread_params = thread_params_vec[exec_id].get();

//     pre_soldering_execs_finished[exec_id] = thread_pool.push([this, thread_params, exec_id, ver_success, inp_from, inp_to, ga_inp_from, ga_inp_to, ga_from, ga_to, &eval_gates_to_blocks, &eval_auths_to_blocks] (int id) {

//       int num_gates = thread_params->num_pre_gates;
//       int num_inputs = thread_params->num_pre_inputs;

//       //Set some often used variables
//       int num_gate_solderings = num_gates * (params.num_bucket - 1) + (num_inputs / 2) * (params.num_inp_bucket - 1);
//       int num_auth_solderings = num_gates * params.num_auth;
//       int num_inp_auth_solderings = num_inputs * (params.num_inp_auth - 1);
//       int num_pre_solderings = 3 * num_gate_solderings + num_auth_solderings + num_inp_auth_solderings;

//       //Receive all preprocessed soldering data and point into this for convenience
//       std::unique_ptr<uint8_t[]> decommited_pre_solderings(std::make_unique<uint8_t[]>(num_pre_solderings * CSEC_BYTES));
//       thread_params->chan.ReceiveBlocking(decommited_pre_solderings.get(), num_pre_solderings * CSEC_BYTES);

//       uint8_t* left_wire_solderings = decommited_pre_solderings.get();
//       uint8_t* right_wire_solderings = left_wire_solderings + CSEC_BYTES * num_gate_solderings;
//       uint8_t* out_wire_solderings = right_wire_solderings + CSEC_BYTES * num_gate_solderings;

//       uint8_t* bucket_auth_solderings = out_wire_solderings + CSEC_BYTES * num_gate_solderings;
//       uint8_t* input_auth_solderings = bucket_auth_solderings + CSEC_BYTES * num_gates * params.num_auth;

//       // Apply the actual solderings and store the indices
//       int curr_head_gate_pos, curr_head_block, curr_head_idx, curr_gate_pos, curr_gate_block, curr_gate_idx, curr_auth_pos, curr_auth_block, curr_auth_idx, curr_head_inp_auth_pos, curr_head_inp_auth_block, curr_head_inp_auth_idx, curr_inp_auth_pos;
//       int solder_gate_pos = 0;
//       int solder_auth_pos = 0;
//       int solder_inp_auth_pos = 0;
//       std::unique_ptr<uint8_t[]> presolder_computed_shares(std::make_unique<uint8_t[]>(num_pre_solderings * CODEWORD_BYTES));

//       //We first loop over all head gates
//       for (int i = ga_from; i < ga_to; ++i) {
//         //Our gate-eval procedure always applies left/right input solderings, so the head gates need to have all-zero left/right solderings. The output soldering is already set, so we ignore this
//         curr_head_gate_pos = i * params.num_bucket;
//         eval_gates_to_blocks.GetExecIDAndIndex(curr_head_gate_pos, curr_head_block, curr_head_idx);

//         //Then all of the gates in this head_gate's bucket (notice j starts at 1).
//         for (int j = 1; j < params.num_bucket; ++j) {
//           curr_gate_pos = curr_head_gate_pos + j;
//           eval_gates_to_blocks.GetExecIDAndIndex(curr_gate_pos, curr_gate_block, curr_gate_idx);

//           //Left soldering
//           std::copy(left_wire_solderings + solder_gate_pos * CSEC_BYTES, left_wire_solderings + solder_gate_pos * CSEC_BYTES + CSEC_BYTES, eval_gates.S_L + curr_gate_pos * CSEC_BYTES);

//           //Decommit shares
//           std::copy(commit_recs[curr_gate_block]->commit_shares[thread_params->left_keys_start + curr_gate_idx], commit_recs[curr_gate_block]->commit_shares[thread_params->left_keys_start + curr_gate_idx] + CODEWORD_BYTES, presolder_computed_shares.get() + solder_gate_pos * CODEWORD_BYTES);
//           XOR_CodeWords(presolder_computed_shares.get() + solder_gate_pos * CODEWORD_BYTES, commit_recs[curr_head_block]->commit_shares[thread_params->left_keys_start + curr_head_idx]);

//           //Right soldering
//           std::copy(right_wire_solderings + solder_gate_pos * CSEC_BYTES, right_wire_solderings + solder_gate_pos * CSEC_BYTES + CSEC_BYTES, eval_gates.S_R + curr_gate_pos * CSEC_BYTES);

//           //Decommit shares
//           std::copy(commit_recs[curr_gate_block]->commit_shares[thread_params->right_keys_start + curr_gate_idx], commit_recs[curr_gate_block]->commit_shares[thread_params->right_keys_start + curr_gate_idx] + CODEWORD_BYTES, presolder_computed_shares.get() + (num_gate_solderings + solder_gate_pos) * CODEWORD_BYTES);
//           XOR_CodeWords(presolder_computed_shares.get() + (num_gate_solderings + solder_gate_pos) * CODEWORD_BYTES, commit_recs[curr_head_block]->commit_shares[thread_params->right_keys_start + curr_head_idx]);

//           //Out soldering
//           XOR_128(eval_gates.S_O + curr_gate_pos * CSEC_BYTES, out_wire_solderings + solder_gate_pos * CSEC_BYTES);

//           //Decommit shares
//           std::copy(commit_recs[curr_gate_block]->commit_shares[thread_params->out_keys_start + curr_gate_idx], commit_recs[curr_gate_block]->commit_shares[thread_params->out_keys_start + curr_gate_idx] + CODEWORD_BYTES, presolder_computed_shares.get() + (2 * num_gate_solderings + solder_gate_pos) * CODEWORD_BYTES);
//           XOR_CodeWords(presolder_computed_shares.get() + (2 * num_gate_solderings + solder_gate_pos) * CODEWORD_BYTES, commit_recs[curr_head_block]->commit_shares[thread_params->out_keys_start + curr_head_idx]);

//           ++solder_gate_pos;
//         }

//         //Then all of the authenticators attached to this head_gate. Here j starts at 0 since all auths are attached to head_gate's output wire
//         for (int j = 0; j < params.num_auth; ++j) {
//           curr_auth_pos = i * params.num_auth + j;
//           eval_auths_to_blocks.GetExecIDAndIndex(curr_auth_pos, curr_auth_block, curr_auth_idx);

//           //Inp_auth soldering
//           std::copy(bucket_auth_solderings + solder_auth_pos * CSEC_BYTES, bucket_auth_solderings + solder_auth_pos * CSEC_BYTES + CSEC_BYTES, eval_auths.S_A + curr_auth_pos * CSEC_BYTES);

//           //Decommit shares
//           std::copy(commit_recs[curr_auth_block]->commit_shares[thread_params->auth_start + curr_auth_idx], commit_recs[curr_auth_block]->commit_shares[thread_params->auth_start + curr_auth_idx] + CODEWORD_BYTES, presolder_computed_shares.get() + (3 * num_gate_solderings + solder_auth_pos) * CODEWORD_BYTES);
//           XOR_CodeWords(presolder_computed_shares.get() + (3 * num_gate_solderings + solder_auth_pos) * CODEWORD_BYTES, commit_recs[curr_head_block]->commit_shares[thread_params->out_keys_start + curr_head_idx]);

//           ++solder_auth_pos;
//         }
//       }

//       //Same as above for gates, but for input buckets
//       for (int i = ga_inp_from; i < ga_inp_to; ++i) {
//         curr_head_gate_pos = params.num_pre_gates * params.num_bucket + i * params.num_inp_bucket;
//         eval_gates_to_blocks.GetExecIDAndIndex(curr_head_gate_pos, curr_head_block, curr_head_idx);

//         //Then all of the gates in this head_gate's bucket (notice j starts at 1).
//         for (int j = 1; j < params.num_inp_bucket; ++j) {
//           curr_gate_pos = curr_head_gate_pos + j;
//           eval_gates_to_blocks.GetExecIDAndIndex(curr_gate_pos, curr_gate_block, curr_gate_idx);

//           //Left soldering
//           std::copy(left_wire_solderings + solder_gate_pos * CSEC_BYTES, left_wire_solderings + solder_gate_pos * CSEC_BYTES + CSEC_BYTES, eval_gates.S_L + curr_gate_pos * CSEC_BYTES);

//           //Decommit shares
//           std::copy(commit_recs[curr_gate_block]->commit_shares[thread_params->left_keys_start + curr_gate_idx], commit_recs[curr_gate_block]->commit_shares[thread_params->left_keys_start + curr_gate_idx] + CODEWORD_BYTES, presolder_computed_shares.get() + solder_gate_pos * CODEWORD_BYTES);
//           XOR_CodeWords(presolder_computed_shares.get() + solder_gate_pos * CODEWORD_BYTES, commit_recs[curr_head_block]->commit_shares[thread_params->left_keys_start + curr_head_idx]);

//           //Right soldering
//           std::copy(right_wire_solderings + solder_gate_pos * CSEC_BYTES, right_wire_solderings + solder_gate_pos * CSEC_BYTES + CSEC_BYTES, eval_gates.S_R + curr_gate_pos * CSEC_BYTES);

//           //Decommit shares
//           std::copy(commit_recs[curr_gate_block]->commit_shares[thread_params->right_keys_start + curr_gate_idx], commit_recs[curr_gate_block]->commit_shares[thread_params->right_keys_start + curr_gate_idx] + CODEWORD_BYTES, presolder_computed_shares.get() + (num_gate_solderings + solder_gate_pos) * CODEWORD_BYTES);
//           XOR_CodeWords(presolder_computed_shares.get() + (num_gate_solderings + solder_gate_pos) * CODEWORD_BYTES, commit_recs[curr_head_block]->commit_shares[thread_params->right_keys_start + curr_head_idx]);

//           //Out soldering
//           XOR_128(eval_gates.S_O + curr_gate_pos * CSEC_BYTES, out_wire_solderings + solder_gate_pos * CSEC_BYTES);

//           //Decommit shares
//           std::copy(commit_recs[curr_gate_block]->commit_shares[thread_params->out_keys_start + curr_gate_idx], commit_recs[curr_gate_block]->commit_shares[thread_params->out_keys_start + curr_gate_idx] + CODEWORD_BYTES, presolder_computed_shares.get() + (2 * num_gate_solderings + solder_gate_pos) * CODEWORD_BYTES);
//           XOR_CodeWords(presolder_computed_shares.get() + (2 * num_gate_solderings + solder_gate_pos) * CODEWORD_BYTES, commit_recs[curr_head_block]->commit_shares[thread_params->out_keys_start + curr_head_idx]);

//           ++solder_gate_pos;
//         }
//       }

//       //Finally we create the indices for input authentication. This is constructed exactly the same as the bucket_solderings as the principle is the same using a head input auth and then solderings onto this all the num_inp_auth-1 other authenticators
//       for (int i = inp_from; i < inp_to; ++i) {
//         curr_head_inp_auth_pos = params.num_pre_gates * params.num_auth + i * params.num_inp_auth;
//         eval_auths_to_blocks.GetExecIDAndIndex(curr_head_inp_auth_pos, curr_head_inp_auth_block, curr_head_inp_auth_idx);

//         for (int j = 1; j < params.num_inp_auth; ++j) {
//           curr_inp_auth_pos = curr_head_inp_auth_pos + j;
//           eval_auths_to_blocks.GetExecIDAndIndex(curr_inp_auth_pos, curr_auth_block, curr_auth_idx);

//           std::copy(input_auth_solderings + solder_inp_auth_pos * CSEC_BYTES, input_auth_solderings + solder_inp_auth_pos * CSEC_BYTES + CSEC_BYTES, eval_auths.S_A + curr_inp_auth_pos * CSEC_BYTES);

//           //Decommit shares
//           std::copy(commit_recs[curr_auth_block]->commit_shares[thread_params->auth_start + curr_auth_idx], commit_recs[curr_auth_block]->commit_shares[thread_params->auth_start + curr_auth_idx] + CODEWORD_BYTES, presolder_computed_shares.get() + (3 * num_gate_solderings + num_auth_solderings + solder_inp_auth_pos) * CODEWORD_BYTES);
//           XOR_CodeWords(presolder_computed_shares.get() + (3 * num_gate_solderings + num_auth_solderings + solder_inp_auth_pos) * CODEWORD_BYTES, commit_recs[curr_head_inp_auth_block]->commit_shares[thread_params->auth_start + curr_head_inp_auth_idx]);

//           ++solder_inp_auth_pos;
//         }
//       }

//       //We end by sending the produced solderings and starting the batch decommit procedure which uses commitments from all executions to build the decommits
//       if (!commit_recs[exec_id]->BatchDecommit(presolder_computed_shares.get(), num_pre_solderings, decommited_pre_solderings.get())) {
//         *ver_success = false;
//         std::cout << "Preprocessed soldering decommit failed" << std::endl;
//       }
    // });
  // }

// //Wait for all executions to finish
//   for (std::future<void>& r : pre_soldering_execs_finished) {
//     r.wait();
//   }
//   auto presolder_end = GET_TIME();

//Check that all executions in both CNC and preprocessed solderings succeeded
  for (std::unique_ptr<bool>& b : thread_ver_successes) {
    if (!*b) {
      throw std::runtime_error("Abort, initial setup failed. Cheating detected");
    }
  }

/////////////////// DEBUG for testing correctness of solderings//////////////
#ifdef DEBUG_SOLDERINGS_INP_BUCKETS
  std::unique_ptr<uint8_t[]> keys_ptr(std::make_unique<uint8_t[]>(CSEC_BYTES * (3 * (params.num_pre_gates + params.num_pre_inputs / 2) + params.num_pre_inputs)));
  uint8_t* keys = keys_ptr.get();
  params.chan.ReceiveBlocking(keys, CSEC_BYTES * (3 * (params.num_pre_gates + params.num_pre_inputs / 2) + params.num_pre_inputs));
  uint8_t out_key[CSEC_BYTES];
  GarblingHandler gh(*thread_params_vec[0]);
  uint8_t* left_key;
  uint8_t* right_key;
  uint8_t* out_key_correct;
  uint8_t* auth_key;

  //Test bucket soldering
  for (int i = 0; i < params.num_pre_gates; ++i) {
    left_key = keys + i * CSEC_BYTES;
    right_key = keys + (params.num_pre_gates + i) * CSEC_BYTES;
    out_key_correct = keys + (2 * params.num_pre_gates + i) * CSEC_BYTES;

    for (int j = 0; j < params.num_bucket; ++j) {
      int curr_gate_pos = i * params.num_bucket + j;
      gh.AllShiftEvaluateGates(eval_gates, curr_gate_pos, left_key, right_key, out_key, eval_gates_ids, 1);

      if (!std::equal(out_key, out_key + CSEC_BYTES, out_key_correct)) {
        std::cout << "gate fail pos:" << i << std::endl;
      }
    }

    for (int j = 0; j < params.num_auth; ++j) {
      int curr_auth_pos = i * params.num_auth + j;
      if (!gh.AllVerifyAuths(eval_auths, curr_auth_pos, out_key, eval_auths_ids, 1)) {
        std::cout << "auth fail pos:" << i << std::endl;
      }
    }
  }

  //Test input buckets soldering
  for (int i = 0; i < params.num_pre_inputs / 2; ++i) {
    left_key = keys + (3 * params.num_pre_gates + i) * CSEC_BYTES;
    right_key = keys + (3 * params.num_pre_gates + params.num_pre_inputs / 2 + i) * CSEC_BYTES;
    out_key_correct = keys + (3 * params.num_pre_gates + 2 * params.num_pre_inputs / 2 + i) * CSEC_BYTES;
    for (int j = 0; j < params.num_inp_bucket; ++j) {
      int curr_gate_pos = params.num_pre_gates * params.num_bucket + i * params.num_inp_bucket + j;
      gh.AllShiftEvaluateGates(eval_gates, curr_gate_pos, left_key, right_key, out_key, eval_gates_ids, 1);
      if (!std::equal(out_key, out_key + CSEC_BYTES, out_key_correct)) {
        std::cout << "input gate fail pos:" << i << std::endl;
      }
    }
  }

  //Test input authentication soldering
  for (int i = 0; i < params.num_pre_inputs; ++i) {
    auth_key = keys + (3 * params.num_pre_gates + 3 * params.num_pre_inputs / 2 + i) * CSEC_BYTES;
    for (int j = 0; j < params.num_inp_auth; ++j) {
      int curr_auth_pos = params.num_pre_gates * params.num_auth + i * params.num_inp_auth + j;
      if (!gh.AllVerifyAuths(eval_auths, curr_auth_pos, auth_key, eval_auths_ids, 1)) {
        std::cout << "inp auth fail pos: " << i << std::endl;
      }
    }
  }
#endif
/////////////////// DEBUG for testing correctness of solderings//////////////

  auto setup_end = GET_TIME();

#ifdef PRINT_COM
  uint64_t bytes_received = params.chan.GetCurrentBytesReceived();
  uint64_t bytes_sent = params.chan.GetCurrentBytesSent();

  params.chan.ResetReceivedBytes();
  params.chan.ResetSentBytes();
  for (std::unique_ptr<Params>& thread_params : thread_params_vec) {
    bytes_received += thread_params->chan.GetCurrentBytesReceived();
    bytes_sent += thread_params->chan.GetCurrentBytesSent();

    thread_params->chan.ResetReceivedBytes();
    thread_params->chan.ResetSentBytes();
  }
  cout << "OT Received " << ot_rec.net.m_vSocket->getRcvCnt() << " bytes" << endl;
  cout << "OT Sent " << ot_rec.net.m_vSocket->getSndCnt() << " bytes" << endl;

  std::cout << "Received " << bytes_received << " Bytes" << std::endl;
  std::cout << "Sent " << bytes_sent << " Bytes" << std::endl;
#endif
  std::vector<std::chrono::duration<long double, std::milli>> durations_res(EVAL_NUM_TIMINGS, std::chrono::duration<long double, std::milli>(0));

  for (int i = 0; i < EVAL_NUM_TIMINGS; ++i) {
    for (int j = 0; j < params.num_execs; ++j)
    {
      durations_res[i] += durations[i][j];
    }
    durations_res[i] = durations_res[i] / params.num_execs;
  }

#ifdef TINY_PRINT
  std::cout << "Avg. Commit: " << durations_res[EVAL_COMMIT_TIME].count() << std::endl;
  std::cout << "Avg. Verleak: " << durations_res[EVAL_VERLEAK_TIME].count() << std::endl;
  std::cout << "Avg. Receive Gates/Auths: " << durations_res[EVAL_RECEIVE_GATES_AUTHS_TIME].count() << std::endl;
  std::cout << "Avg. CnC: " << durations_res[EVAL_CNC_TIME].count() << std::endl;
  PRINT_TIME(presolder_end, presolder_begin, "PRE_SOLDER");
  PRINT_TIME(setup_end, setup_begin, "SETUP_TOTAL");
#endif
}

void TinyEvaluator::Offline(std::vector<Circuit*>& circuits, int top_num_execs) {
//   // int top_num_execs = std::min((int)circuits.size(), params.num_execs);
//   std::vector<std::future<void>> top_soldering_execs_finished(top_num_execs);
//   std::vector<std::unique_ptr<bool>> thread_ver_successes;

//   //Split the number of preprocessed gates and inputs into top_num_execs executions
//   std::vector<int> circuits_from, circuits_to;
//   PartitionBufferFixedNum(circuits_from, circuits_to, top_num_execs, circuits.size());

//   int num_gates_needed = 0;
//   int num_inp_gates_needed = 0;
//   int num_inps_needed = 0;
//   int num_outs_needed = 0;
//   for (int i = 0; i < circuits.size(); ++i) {
//     gates_offset.emplace_back(num_gates_used + num_gates_needed);
//     inp_gates_offset.emplace_back(num_inputs_used / 2 + num_inp_gates_needed);
//     inputs_offset.emplace_back(num_inputs_used + num_inps_needed);
//     outputs_offset.emplace_back(num_outputs_used + num_outs_needed);
//     num_gates_needed += circuits[i]->num_and_gates;
//     num_inp_gates_needed += circuits[i]->num_const_inp_wires / 2;
//     num_inps_needed += circuits[i]->num_inp_wires;
//     num_outs_needed += circuits[i]->num_out_wires;
//   }

//   if ((params.num_pre_gates - num_gates_used) < num_gates_needed) {
//     throw std::runtime_error("Not enough garbled gates");
//   } else {
//     num_gates_used += num_gates_needed;
//   }

//   //Due to the way we choose our parameters, if there are enough num_inps, then there are also enough for inp_gates.
//   if ((params.num_pre_inputs - num_inputs_used) < num_inps_needed) {
//     throw std::runtime_error("Not enough input authenticators");
//   } else {
//     num_inputs_used += num_inps_needed;
//   }

//   if ((params.num_pre_outputs - num_outputs_used) < num_outs_needed) {
//     throw std::runtime_error("Not enough output wires");
//   } else {
//     num_outputs_used += num_outs_needed;
//   }

//   //Setup maps from eval_gates and eval_auths to commit_block and inner block commit index. Needed to construct decommits that span all executions
//   IDMap eval_gates_to_blocks(eval_gates_ids, thread_params_vec[0]->Q + thread_params_vec[0]->A, thread_params_vec[0]->out_keys_start);
//   IDMap eval_auths_to_blocks(eval_auths_ids, thread_params_vec[0]->Q + thread_params_vec[0]->A, thread_params_vec[0]->auth_start);

//   auto topo_soldering_begin = GET_TIME();
//   for (int exec_id = 0; exec_id < top_num_execs; ++exec_id) {
//     int circ_from = circuits_from[exec_id];
//     int circ_to = circuits_to[exec_id];
//     Params* thread_params = thread_params_vec[exec_id].get();
//     thread_ver_successes.emplace_back(std::make_unique<bool>(true));
//     bool* ver_success = thread_ver_successes[exec_id].get();

//     top_soldering_execs_finished[exec_id] = thread_pool.push([this, thread_params, exec_id, circ_from, circ_to, &circuits, &eval_gates_to_blocks, &eval_auths_to_blocks, ver_success] (int id) {

//       for (int c = circ_from; c < circ_to; ++c) {
//         Circuit* circuit = circuits[c];
//         int gate_offset = gates_offset[c];
//         int inp_gate_offset = inp_gates_offset[c];
//         int inp_offset = inputs_offset[c];
//         int num_top_solderings = 2 * circuit->num_and_gates + circuit->num_const_inp_wires; //the 2* factor cancels out as we can check two inputs pr. input bucket.

//         std::unique_ptr<uint8_t[]> topsolder_computed_shares(std::make_unique<uint8_t[]>(num_top_solderings * CODEWORD_BYTES + circuit->num_wires * CODEWORD_BYTES));
//         uint8_t* topsolder_computed_shares_tmp = topsolder_computed_shares.get() + num_top_solderings * CODEWORD_BYTES;

//         int curr_auth_inp_head_pos, curr_inp_head_block, curr_inp_head_idx, curr_head_pos, curr_head_block, curr_head_idx;
//         for (int i = 0; i < circuit->num_inp_wires; ++i) {
//           curr_auth_inp_head_pos = params.num_pre_gates * thread_params->num_auth + (inp_offset + i) * thread_params->num_inp_auth;
//           eval_auths_to_blocks.GetExecIDAndIndex(curr_auth_inp_head_pos, curr_inp_head_block, curr_inp_head_idx);

//           //Build decommit_info
//           std::copy(commit_recs[curr_inp_head_block]->commit_shares[thread_params->auth_start + curr_inp_head_idx], commit_recs[curr_inp_head_block]->commit_shares[thread_params->auth_start + curr_inp_head_idx] + CODEWORD_BYTES, topsolder_computed_shares_tmp + i * CODEWORD_BYTES);
//         }

//         int left_inp_start = circuit->num_and_gates;
//         int right_inp_start = 2 * circuit->num_and_gates + circuit->num_const_inp_wires / 2;
//         for (int i = 0; i < circuit->num_const_inp_wires / 2; ++i) {
//           curr_head_pos = params.num_pre_gates * params.num_bucket + (inp_gate_offset + i) * params.num_inp_bucket;
//           eval_gates_to_blocks.GetExecIDAndIndex(curr_head_pos, curr_head_block, curr_head_idx);

//           std::copy(commit_recs[curr_head_block]->commit_shares[thread_params->left_keys_start + curr_head_idx], commit_recs[curr_head_block]->commit_shares[thread_params->left_keys_start + curr_head_idx] + CODEWORD_BYTES, topsolder_computed_shares.get() + (left_inp_start + i) * CODEWORD_BYTES);
//           XOR_CodeWords(topsolder_computed_shares.get() + (left_inp_start + i) * CODEWORD_BYTES, topsolder_computed_shares_tmp + i * CODEWORD_BYTES);
//           std::copy(commit_recs[curr_head_block]->commit_shares[thread_params->right_keys_start + curr_head_idx], commit_recs[curr_head_block]->commit_shares[thread_params->right_keys_start + curr_head_idx] + CODEWORD_BYTES, topsolder_computed_shares.get() + (right_inp_start + i) * CODEWORD_BYTES);
//           XOR_CodeWords(topsolder_computed_shares.get() + (right_inp_start + i) * CODEWORD_BYTES, topsolder_computed_shares_tmp + (circuit->num_const_inp_wires / 2 + i) * CODEWORD_BYTES);
//         }

//         int curr_and_gate = 0;
//         int left_gate_start = 0;
//         int right_gate_start = circuit->num_and_gates + circuit->num_const_inp_wires / 2;
//         for (int i = 0; i < circuit->num_gates; ++i) {
//           Gate g = circuit->gates[i];
//           if (g.type == NOT) {
//             //Build decommit_info
//             std::copy(topsolder_computed_shares_tmp + g.left_wire * CODEWORD_BYTES, topsolder_computed_shares_tmp + g.left_wire * CODEWORD_BYTES + CODEWORD_BYTES, topsolder_computed_shares_tmp + g.out_wire * CODEWORD_BYTES);

//             XOR_CodeWords(topsolder_computed_shares_tmp + g.out_wire * CODEWORD_BYTES, commit_recs[0]->commit_shares[thread_params->delta_pos]);

//           } else if (g.type == XOR) {
//             //Build decommit_info
//             std::copy(topsolder_computed_shares_tmp + g.left_wire * CODEWORD_BYTES, topsolder_computed_shares_tmp + g.left_wire * CODEWORD_BYTES + CODEWORD_BYTES, topsolder_computed_shares_tmp + g.out_wire * CODEWORD_BYTES);

//             XOR_CodeWords(topsolder_computed_shares_tmp + g.out_wire * CODEWORD_BYTES, topsolder_computed_shares_tmp + g.right_wire * CODEWORD_BYTES);

//           } else if (g.type == AND) {
//             curr_head_pos = (gate_offset + curr_and_gate) * thread_params->num_bucket;
//             eval_gates_to_blocks.GetExecIDAndIndex(curr_head_pos, curr_head_block, curr_head_idx);
//             //Build decommit_info
//             std::copy(commit_recs[curr_head_block]->commit_shares[thread_params->out_keys_start + curr_head_idx], commit_recs[curr_head_block]->commit_shares[thread_params->out_keys_start + curr_head_idx] + CODEWORD_BYTES, topsolder_computed_shares_tmp + g.out_wire * CODEWORD_BYTES);

//             std::copy(commit_recs[curr_head_block]->commit_shares[thread_params->left_keys_start + curr_head_idx], commit_recs[curr_head_block]->commit_shares[thread_params->left_keys_start + curr_head_idx] + CODEWORD_BYTES, topsolder_computed_shares.get() + (left_gate_start + curr_and_gate) * CODEWORD_BYTES);

//             std::copy(commit_recs[curr_head_block]->commit_shares[thread_params->right_keys_start + curr_head_idx], commit_recs[curr_head_block]->commit_shares[thread_params->right_keys_start + curr_head_idx] + CODEWORD_BYTES, topsolder_computed_shares.get() + (right_gate_start + curr_and_gate) * CODEWORD_BYTES);

//             XOR_CodeWords(topsolder_computed_shares.get() + (left_gate_start + curr_and_gate) * CODEWORD_BYTES, topsolder_computed_shares_tmp + g.left_wire * CODEWORD_BYTES);

//             XOR_CodeWords(topsolder_computed_shares.get() + (right_gate_start + curr_and_gate) * CODEWORD_BYTES, topsolder_computed_shares_tmp + g.right_wire * CODEWORD_BYTES);

//             ++curr_and_gate;
//           }
//         }

//         std::unique_ptr<uint8_t[]> topological_solderings(std::make_unique<uint8_t[]>(num_top_solderings * CSEC_BYTES));

//         thread_params->chan.ReceiveBlocking(topological_solderings.get(), num_top_solderings * CSEC_BYTES);
//         if (!commit_recs[exec_id]->BatchDecommit(topsolder_computed_shares.get(), num_top_solderings, topological_solderings.get())) {
//           *ver_success = false;
//           std::cout << "Topological soldering decommit failed" << std::endl;
//         }

//         for (int i = 0; i < circuit->num_and_gates; ++i) {
//           for (int j = 0; j < thread_params->num_bucket; ++j) {
//             curr_head_pos = (gate_offset + i) * thread_params->num_bucket + j;
//             XOR_128(eval_gates.S_L + curr_head_pos * CSEC_BYTES, topological_solderings.get() + (left_gate_start + i) * CSEC_BYTES);

//             XOR_128(eval_gates.S_R + curr_head_pos * CSEC_BYTES, topological_solderings.get() + (right_gate_start + i) * CSEC_BYTES);
//           }
//         }

//         for (int i = 0; i < circuit->num_const_inp_wires / 2; ++i) {
//           for (int j = 0; j < thread_params->num_inp_bucket; ++j) {
//             curr_head_pos = params.num_pre_gates * params.num_bucket + (inp_gate_offset + i) * params.num_inp_bucket;
//             XOR_128(eval_gates.S_L + (curr_head_pos + j) * CSEC_BYTES, topological_solderings.get() + (left_inp_start + i) * CSEC_BYTES);

//             XOR_128(eval_gates.S_R + (curr_head_pos + j) * CSEC_BYTES, topological_solderings.get() + (right_inp_start + i) * CSEC_BYTES);
//           }
//         }
//       }
//     });
//   }

//   for (std::future<void>& r : top_soldering_execs_finished) {
//     r.wait();
//   }
//   auto topo_soldering_end = GET_TIME();

//   for (std::unique_ptr<bool>& b : thread_ver_successes) {
//     if (!*b) {
//       throw std::runtime_error("Abort, topological soldering failed. Cheating detected");
//     }
//   }

// #ifdef PRINT_COM
//   uint64_t bytes_received = params.chan.GetCurrentBytesReceived();
//   uint64_t bytes_sent = params.chan.GetCurrentBytesSent();

//   params.chan.ResetReceivedBytes();
//   params.chan.ResetSentBytes();
//   for (std::unique_ptr<Params>& thread_params : thread_params_vec) {
//     bytes_received += thread_params->chan.GetCurrentBytesReceived();
//     bytes_sent += thread_params->chan.GetCurrentBytesSent();

//     thread_params->chan.ResetReceivedBytes();
//     thread_params->chan.ResetSentBytes();
//   }
//   std::cout << "Received " << bytes_received << " Bytes" << std::endl;
//   std::cout << "Sent " << bytes_sent << " Bytes" << std::endl;
// #endif

// #ifdef TINY_PRINT
//   PRINT_TIME(topo_soldering_end, topo_soldering_begin, "TOP_SOLDER");
// #endif
}

void TinyEvaluator::Online(std::vector<Circuit*>& circuits, std::vector<uint8_t*>& inputs, std::vector<uint8_t*>& outputs, int eval_num_execs) {

//   std::vector<std::future<void>> online_execs_finished(eval_num_execs);
//   std::vector<int> circuits_from, circuits_to;

//   PartitionBufferFixedNum(circuits_from, circuits_to, eval_num_execs, circuits.size());

//   IDMap eval_auths_to_blocks(eval_auths_ids, thread_params_vec[0]->Q + thread_params_vec[0]->A, thread_params_vec[0]->auth_start);
//   IDMap eval_gates_to_blocks(eval_gates_ids, thread_params_vec[0]->Q + thread_params_vec[0]->A, thread_params_vec[0]->out_keys_start);

//   for (int exec_id = 0; exec_id < eval_num_execs; ++exec_id) {

//     int circ_from = circuits_from[exec_id];
//     int circ_to = circuits_to[exec_id];
//     Params* thread_params = thread_params_vec[exec_id].get();

//     online_execs_finished[exec_id] = thread_pool.push([this, thread_params, exec_id, circ_from, circ_to, &circuits, &inputs, &outputs, &eval_gates_to_blocks, &eval_auths_to_blocks] (int id) {

//       Circuit* circuit;
//       Gate g;
//       uint8_t* eval_input;
//       uint8_t* eval_outputs;
//       uint8_t* const_inp_keys;
//       uint8_t* eval_inp_keys;
//       uint8_t* out_decommit_values;
//       uint8_t* eval_computed_shares_inp;
//       uint8_t* eval_computed_shares_out;
//       uint8_t* decommit_shares_inp_0;
//       uint8_t* decommit_shares_inp_1;
//       uint8_t* decommit_shares_out_0;
//       uint8_t* decommit_shares_out_1;
//       uint8_t* e;
//       __m128i intrin_outs[thread_params->num_bucket];
//       __m128i intrin_auths[thread_params->num_auth];
//       int bucket_score[thread_params->num_bucket];
//       bool all_equal;

//       int gate_offset, inp_gate_offset, inp_offset, out_offset, curr_and_gate;
//       int curr_auth_inp_head_pos, curr_inp_head_block, curr_inp_head_idx, curr_head_pos, curr_head_block, curr_head_idx, curr_output_pos, curr_output_block, curr_output_idx;

//       GarblingHandler gh(*thread_params);
//       int curr_input, curr_output, ot_commit_block, commit_id, chosen_val_id;
//       for (int c = circ_from; c < circ_to; ++c) {
//         auto t_0 = GET_TIME();
//         circuit = circuits[c];
//         eval_input = inputs[c];
//         eval_outputs = outputs[c];
//         gate_offset = gates_offset[c];
//         inp_gate_offset = inp_gates_offset[c];
//         inp_offset = inputs_offset[c];
//         out_offset = outputs_offset[c];

//         uint8_t* ot_input = new uint8_t[2 * BITS_TO_BYTES(circuit->num_eval_inp_wires) + (circuit->num_eval_inp_wires + circuit->num_out_wires) * CODEWORD_BYTES + circuit->num_const_inp_wires * CSEC_BYTES + (circuit->num_eval_inp_wires + circuit->num_out_wires) * (CODEWORD_BYTES + 2 * CSEC_BYTES)];

//         e = ot_input + BITS_TO_BYTES(circuit->num_eval_inp_wires);

//         eval_computed_shares_inp = e + BITS_TO_BYTES(circuit->num_eval_inp_wires);
//         eval_computed_shares_out = eval_computed_shares_inp + circuit->num_eval_inp_wires * CODEWORD_BYTES;

//         const_inp_keys = eval_computed_shares_out + circuit->num_out_wires * CODEWORD_BYTES;

//         decommit_shares_inp_0 = const_inp_keys + circuit->num_const_inp_wires * CSEC_BYTES;
//         decommit_shares_inp_1 = decommit_shares_inp_0 + circuit->num_eval_inp_wires * CODEWORD_BYTES;

//         eval_inp_keys = decommit_shares_inp_1 + circuit->num_eval_inp_wires * CSEC_BYTES;

//         decommit_shares_out_0 = eval_inp_keys + circuit->num_eval_inp_wires * CSEC_BYTES;
//         decommit_shares_out_1 = decommit_shares_out_0 + circuit->num_out_wires * CODEWORD_BYTES;

//         out_decommit_values = decommit_shares_out_1 + circuit->num_eval_inp_wires * CSEC_BYTES;

//         uint32_t num_receiving_bytes_inp = circuit->num_const_inp_wires * CSEC_BYTES + circuit->num_eval_inp_wires * (CODEWORD_BYTES + CSEC_BYTES);
//         uint32_t num_receiving_bytes_out = circuit->num_out_wires * (CODEWORD_BYTES + CSEC_BYTES);

//         for (int i = 0; i < circuit->num_eval_inp_wires; ++i) {
//           curr_input = (inp_offset + i);
//           if (GetBit(curr_input, ot_rec.choices_outer.get())) {
//             SetBit(i, 1, ot_input);
//           } else {
//             SetBit(i, 0, ot_input);
//           }
//         }

//         XOR_UINT8_T(e, eval_input, ot_input, BITS_TO_BYTES(circuit->num_eval_inp_wires));
//         thread_params->chan.Send(e, BITS_TO_BYTES(circuit->num_eval_inp_wires));

//         auto t0 = GET_TIME();

//         __m128i* intrin_values = new __m128i[circuit->num_wires]; //using raw pointer due to ~25% increase in overall performance. Since the online phase is so computationally efficient even the slightest performance hit is immediately seen. It does not matter in the others phases as they operation on a very different running time scale.

//         for (int i = 0; i < circuit->num_eval_inp_wires; ++i) {
//           curr_input = (inp_offset + i);
//           ot_commit_block = curr_input / thread_params->num_pre_inputs;
//           commit_id = thread_params->ot_chosen_start + curr_input % thread_params->num_pre_inputs;

//           std::copy(commit_recs[ot_commit_block]->commit_shares[commit_id], commit_recs[ot_commit_block]->commit_shares[commit_id] + CODEWORD_BYTES, eval_computed_shares_inp + i * CODEWORD_BYTES);

//           //Add the input key
//           curr_auth_inp_head_pos = params.num_pre_gates * thread_params->num_auth + (inp_offset + circuit->num_const_inp_wires + i) * thread_params->num_inp_auth;
//           eval_auths_to_blocks.GetExecIDAndIndex(curr_auth_inp_head_pos, curr_inp_head_block, curr_inp_head_idx);
//           XOR_CodeWords(eval_computed_shares_inp + i * CODEWORD_BYTES, commit_recs[curr_inp_head_block]->commit_shares[thread_params->auth_start + curr_inp_head_idx]);

//           if (GetBit(i, e)) {
//             XOR_CodeWords(eval_computed_shares_inp + i * CODEWORD_BYTES, commit_recs[ot_commit_block]->commit_shares[thread_params->delta_pos]);
//           }
//         }

//         auto t = GET_TIME();
//         thread_params->chan.ReceiveBlocking(const_inp_keys, num_receiving_bytes_inp);
//         auto t2 = GET_TIME();

//         decommit_shares_inp_0 = const_inp_keys + circuit->num_const_inp_wires * CSEC_BYTES;
//         decommit_shares_inp_1 = decommit_shares_inp_0 + circuit->num_eval_inp_wires * CODEWORD_BYTES;
//         decommit_shares_out_0 = decommit_shares_inp_1 + circuit->num_eval_inp_wires * CSEC_BYTES;
//         decommit_shares_out_1 = decommit_shares_out_0 + circuit->num_out_wires * CODEWORD_BYTES;
//         if (!VerifyDecommits(decommit_shares_inp_0, decommit_shares_inp_1, eval_computed_shares_inp, eval_inp_keys, rot_choices.get(), commit_recs[exec_id]->code.get(), circuit->num_eval_inp_wires)) {
//           throw std::runtime_error("Abort: Wrong eval keys sent!");
//         }

//         auto t3 = GET_TIME();

//         for (int i = 0; i < circuit->num_eval_inp_wires; ++i) {
//           curr_input = (inp_offset + i);
//           ot_commit_block = curr_input / thread_params->num_pre_inputs;
//           chosen_val_id = curr_input % thread_params->num_pre_inputs;
//           XOR_128(eval_inp_keys + i * CSEC_BYTES, commit_recs[ot_commit_block]->chosen_commit_values.get() + chosen_val_id * CSEC_BYTES);

//           //XOR out lsb(K^i_0)
//           uint8_t lsb_zero_key = GetLSB(eval_inp_keys + i * CSEC_BYTES) ^ GetBit(params.num_pre_outputs + inp_offset + i, verleak_bits.get()) ^ GetBit(i, e);

//           //XOR out K^i_y_i
//           XOR_128(eval_inp_keys + i * CSEC_BYTES, ot_rec.response_outer.get() + curr_input * CSEC_BYTES);

//           //Check using lsb(K^i_0) that we received the correct key according to eval_input
//           if ((GetLSB(eval_inp_keys + i * CSEC_BYTES) ^ lsb_zero_key) != GetBit(i, eval_input)) {
//             throw std::runtime_error("Abort: Wrong eval value keys sent!");
//           }

//           intrin_values[circuit->num_const_inp_wires + i] = _mm_lddqu_si128((__m128i *) (eval_inp_keys + i * CSEC_BYTES));
//         }

//         auto t4 = GET_TIME();

//         //Ensure that constructor sends valid keys. Implemented different than in paper as we here require that ALL input authenticators accept. This has no influence on security as if we abort here it does not leak anything about the evaluators input. Also, the sender knows if the evaluator is going to abort before sending the bad keys, so it leaks nothing.
//         for (int i = 0; i < circuit->num_const_inp_wires; ++i) {
//           intrin_values[i] = _mm_lddqu_si128((__m128i *) (const_inp_keys + i * CSEC_BYTES));
//         }

//         for (int i = 0; i < circuit->num_inp_wires; ++i) {
//           curr_auth_inp_head_pos = params.num_pre_gates * thread_params->num_auth + (inp_offset + i) * thread_params->num_inp_auth;

//           for (int j = 0; j < thread_params->num_inp_auth; ++j) {
//             if (!IntrinVerifyAuths(eval_auths, curr_auth_inp_head_pos + j, intrin_values[i], eval_auths_ids[curr_auth_inp_head_pos + j], gh.key_schedule)) {
//               throw std::runtime_error("Abort: Inp auth fail!");
//             }
//           }
//         }

// /////////////////////////////// DEBUG Input buckets////////////////////////////
// #ifdef DEBUG_SOLDERINGS_INP_BUCKETS
//         __m128i out_keys[2];
//         for (int i = 0; i < circuit->num_const_inp_wires / 2; ++i) {
//           int curr_inp_gate_head = params.num_pre_gates * params.num_bucket + (inp_gate_offset + i) * params.num_inp_bucket;
//           IntrinShiftEvaluateGates(eval_gates, curr_inp_gate_head, intrin_values[i], intrin_values[circuit->num_const_inp_wires / 2 + i], out_keys[0], eval_gates_ids[curr_inp_gate_head], gh.key_schedule);

//           for (int j = 1; j < params.num_inp_bucket; ++j) {
//             IntrinShiftEvaluateGates(eval_gates, curr_inp_gate_head + j, intrin_values[i], intrin_values[circuit->num_const_inp_wires / 2 + i], out_keys[1], eval_gates_ids[curr_inp_gate_head + j], gh.key_schedule);
//             if (!compare128(out_keys[0], out_keys[1])) {
//               std::cout << "input gate fail pos:" << i << std::endl;
//             }
//           }
//         }
// #endif
// /////////////////////////////// DEBUG Input buckets////////////////////////////

//         auto t5 = GET_TIME();
//         curr_and_gate = 0;
//         for (int i = 0; i < circuit->num_gates; ++i) {
//           g = circuit->gates[i];
//           if (g.type == NOT) {
//             intrin_values[g.out_wire] = intrin_values[g.left_wire];
//           } else if (g.type == XOR) {
//             intrin_values[g.out_wire] = _mm_xor_si128(intrin_values[g.left_wire], intrin_values[g.right_wire]);
//           } else if (g.type == AND) {
//             curr_head_pos = (gate_offset + curr_and_gate) * thread_params->num_bucket;

//             IntrinShiftEvaluateGates(eval_gates, curr_head_pos, intrin_values[g.left_wire], intrin_values[g.right_wire], intrin_values[g.out_wire], eval_gates_ids[curr_head_pos], gh.key_schedule);

//             all_equal = true;
//             for (int j = 1; j < thread_params->num_bucket; ++j) {
//               IntrinShiftEvaluateGates(eval_gates, curr_head_pos + j, intrin_values[g.left_wire], intrin_values[g.right_wire], intrin_outs[j - 1], eval_gates_ids[curr_head_pos + j], gh.key_schedule);
//               all_equal &= compare128(intrin_values[g.out_wire], intrin_outs[j - 1]);
//             }

//             if (!all_equal) {
//               std::cout << "all outputs not equal for " << curr_and_gate << std::endl;
//               std::fill(bucket_score, bucket_score + thread_params->num_bucket * sizeof(uint32_t), 0);
//               intrin_outs[thread_params->num_bucket - 1] = intrin_values[g.out_wire]; //now all keys are in intrin_outs.
//               intrin_auths[0] = intrin_outs[0];
//               ++bucket_score[0];
//               int candidates = 1;
//               for (int j = 1; j < thread_params->num_bucket; ++j) {
//                 int comp = 0;
//                 for (int k = 0; k < candidates; k++) {
//                   comp = !compare128(intrin_outs[j], intrin_outs[k]);
//                   if (comp == 0) {
//                     ++bucket_score[k];
//                     break;
//                   }
//                 }
//                 if (comp != 0) {
//                   intrin_auths[candidates] = intrin_outs[j];
//                   ++candidates;
//                 }
//               }

//               //Check the candidates
//               for (int j = 0; j < thread_params->num_auth; j++) {
//                 curr_auth_inp_head_pos = params.num_pre_gates * thread_params->num_auth + (inp_offset + curr_and_gate) * thread_params->num_inp_auth;

//                 for (uint32_t k = 0; k < candidates; k++) {
//                   int res = IntrinVerifyAuths(eval_auths, curr_auth_inp_head_pos + k, intrin_auths[k], eval_auths_ids[curr_auth_inp_head_pos + k], gh.key_schedule);
//                   if (res == 1) { // The key is good
//                     ++bucket_score[k];
//                   }
//                 }
//               }
//               // Find the winner
//               int winner_idx = -1;
//               int curr_high_score = -1;
//               for (int j = 0; j < candidates; j++) {
//                 if (bucket_score[j] > curr_high_score) {
//                   winner_idx = j;
//                   curr_high_score = bucket_score[j];
//                 }
//               }

//               intrin_values[g.out_wire] = intrin_auths[winner_idx];
//             }
//             ++curr_and_gate;
//           }
//         }

//         for (int i = 0; i < circuit->num_out_wires; ++i) {
//           curr_output = (out_offset + i);
//           ot_commit_block = curr_output / thread_params->num_pre_outputs;
//           commit_id = thread_params->out_lsb_blind_start + curr_output % thread_params->num_pre_outputs;

//           std::copy(commit_recs[ot_commit_block]->commit_shares[commit_id], commit_recs[ot_commit_block]->commit_shares[commit_id] + CODEWORD_BYTES, eval_computed_shares_out + i * CODEWORD_BYTES);

//           //Add the output key
//           curr_output_pos = (gate_offset + circuit->num_and_gates - circuit->num_out_wires + i) * thread_params->num_bucket;
//           eval_gates_to_blocks.GetExecIDAndIndex(curr_output_pos, curr_output_block, curr_output_idx);

//           XOR_CodeWords(eval_computed_shares_out + i * CODEWORD_BYTES, commit_recs[curr_output_block]->commit_shares[thread_params->out_keys_start + curr_output_idx]);
//         }

//         thread_params->chan.ReceiveBlocking(decommit_shares_out_0, num_receiving_bytes_out);

//         if (!VerifyDecommits(decommit_shares_out_0, decommit_shares_out_1, eval_computed_shares_out, out_decommit_values, rot_choices.get(), commit_recs[exec_id]->code.get(), circuit->num_out_wires)) {
//           throw std::runtime_error("Abort: Wrong eval keys sent!");
//         }

//         auto t6 = GET_TIME();
//         for (int i = 0; i < circuit->num_out_wires; ++i) {
//           SetBit(i, GetLSB(out_decommit_values + i * CSEC_BYTES) ^ GetBit(out_offset + i, verleak_bits.get()) ^ GetLSB(intrin_values[circuit->num_wires - circuit->num_out_wires + i]), eval_outputs);
//         }

//         auto t7 = GET_TIME();

//         delete[] ot_input;
//         delete[] intrin_values;

// #ifdef TINY_PRINT
//         //Could also report average as in preprocessing
//         if ((exec_id == 0) && (c == 0)) {
//           PRINT_TIME_NANO(t0, t_0, "inp_prep");
//           PRINT_TIME_NANO(t, t0, "commit_share");
//           PRINT_TIME_NANO(t2, t, "key_wait");
//           PRINT_TIME_NANO(t3, t2, "commit");
//           PRINT_TIME_NANO(t4, t3, "eval inp");
//           PRINT_TIME_NANO(t5, t4, "const inp");
//           PRINT_TIME_NANO(t6, t5, "eval circ");
//           PRINT_TIME_NANO(t7, t6, "output decoding");
//         }
// #endif
//       }
//     });
//   }

//   for (std::future<void>& r : online_execs_finished) {
//     r.wait();
//   }

// #ifdef PRINT_COM
//   uint64_t bytes_received = params.chan.GetCurrentBytesReceived();
//   uint64_t bytes_sent = params.chan.GetCurrentBytesSent();

//   params.chan.ResetReceivedBytes();
//   params.chan.ResetSentBytes();
//   for (std::unique_ptr<Params>& thread_params : thread_params_vec) {
//     bytes_received += thread_params->chan.GetCurrentBytesReceived();
//     bytes_sent += thread_params->chan.GetCurrentBytesSent();

//     thread_params->chan.ResetReceivedBytes();
//     thread_params->chan.ResetSentBytes();
//   }
//   std::cout << "Received " << bytes_received << " Bytes" << std::endl;
//   std::cout << "Sent " << bytes_sent << " Bytes" << std::endl;
// #endif
}

//The below function is essentially a mix of the three CommitRec member functions ConsistencyCheck, BatchDecommit and VerifyTransposedDecommits.
bool TinyEvaluator::BatchDecommitLSB(CommitReceiver* commit_rec, uint8_t decommit_shares[], int num_values, uint8_t values[]) {

  // uint8_t ver_leak_challenge[BITS_TO_BYTES(commit_rec->col_dim_single) + CSEC_BYTES];
  // commit_rec->params.rnd.GenRnd(ver_leak_challenge, 2 * CSEC_BYTES);
  // commit_rec->params.chan.Send(ver_leak_challenge, 2 * CSEC_BYTES);

  // uint8_t* delta_chal = ver_leak_challenge;
  // uint8_t* alpha_seed = ver_leak_challenge + CSEC_BYTES;

  // //Preprocess DELTA
  // int delta_matrix_size = BITS_TO_BYTES(commit_rec->row_dim * commit_rec->col_dim_single); //delta_matrix_size is 8x smaller than transpose_matrix_size.

  // std::unique_ptr<uint8_t[]> delta_matrix_tmp(std::make_unique<uint8_t[]>(2 * delta_matrix_size));

  // //Load Delta share into the matrix to be transposed. However only if the bit is set in delta_chal. The remaining columns are left as 0. This reflects adding Delta or not.
  // for (int i = 0; i < commit_rec->col_dim_single; ++i) {
  //   if (GetBit(i, delta_chal)) {
  //     std::copy(commit_rec->commit_shares[commit_rec->params.delta_pos], commit_rec->commit_shares[commit_rec->params.delta_pos] + CODEWORD_BYTES, delta_matrix_tmp.get() + delta_matrix_size + i * commit_rec->row_dim_bytes);
  //   }
  // }

  // //Transpose the blocks containing the Delta share
  // transpose_128_320(delta_matrix_tmp.get() + delta_matrix_size, delta_matrix_tmp.get(), 1);


  // //The below proceeds more in line with normal batch decommit/consistency check.

  // //Setup all registers for calculation the linear combinations. Will end up with SSEC linear combinations. The final_values_result will hold the linear combinations of the postulated values.
  // uint8_t final_result[CODEWORD_BYTES * SSEC];
  // uint8_t final_values_result[CSEC_BYTES * SSEC];

  // //res_tmp is twice as large as we do not do degree reduction until the very end, so we need to accumulate a larger intermediate value.
  // __m128i res_tmp[2][CODEWORD_BITS];
  // __m128i res_values_tmp[2][AES_BITS];
  // __m128i res_total[CODEWORD_BITS];
  // for (int i = 0; i < CODEWORD_BITS; ++i) {
  //   res_tmp[0][i] = _mm_setzero_si128();
  //   res_tmp[1][i] = _mm_setzero_si128();
  //   res_total[i] = _mm_lddqu_si128((__m128i*) (delta_matrix_tmp.get() + i * AES_BYTES)); //Apply DELTA shares
  //   if (i < AES_BITS) {
  //     res_values_tmp[0][i] = _mm_setzero_si128();
  //     res_values_tmp[1][i] = _mm_setzero_si128();
  //   }
  // }

  // __m128i val;
  // __m128i val_result[2];

  // //Used for transposing the commitment shares
  // std::unique_ptr<uint8_t[]> matrices_tmp(std::make_unique<uint8_t[]>(2 * commit_rec->transpose_matrix_size));

  // //Used for transposing the values
  // CBitVector matrix_buffer_values;
  // matrix_buffer_values.AttachBuf((uint8_t*) malloc(commit_rec->transpose_matrix_values_size), commit_rec->transpose_matrix_values_size);

  // //Load the initial challenge element
  // __m128i alpha = _mm_lddqu_si128((__m128i *) alpha_seed);

  // //Compute number of check_blocks needed in total for num_values
  // int num_check_blocks = CEIL_DIVIDE(num_values, commit_rec->col_dim);

  // //For each check_block we load the shares in column-major order and then transpose to get to row-major order so we can address AES_BITS values entry-wise at a time.
  // for (int j = 0; j < num_check_blocks; ++j) {
  //   //As we only compare bits we set the vector to all 0 and load the bits one by one. Therefore need to reset it each time.
  //   std::fill(matrix_buffer_values.GetArr(), matrix_buffer_values.GetArr() + commit_rec->transpose_matrix_values_size, 0);

  //   //Load block
  //   for (int i = 0; i < commit_rec->col_dim; ++i) {

  //     //First we check if we are done, ie we are in the last block and it is not entirely filled up.
  //     int num_check_index = j * commit_rec->col_dim + i;
  //     if (num_check_index < num_values) {
  //       //We copy the correct share into our current block matrix
  //       std::copy(decommit_shares + num_check_index * CODEWORD_BYTES, decommit_shares + num_check_index * CODEWORD_BYTES + CODEWORD_BYTES, matrices_tmp.get() + commit_rec->transpose_matrix_size + i * commit_rec->row_dim_bytes);

  //       //We copy the current bit value into our current block matrix
  //       SetBit(AES_BITS - 1, GetBit(num_check_index, values), matrix_buffer_values.GetArr() + i * commit_rec->row_dim_values_bytes);

  //     } else {
  //       //This pads the last block with 0 rows. Can be optimized by only doing this once. However need to calculate the remaining number of "rows" to zero out then.
  //       std::fill(matrices_tmp.get() + commit_rec->transpose_matrix_size + i * commit_rec->row_dim_bytes, matrices_tmp.get() + commit_rec->transpose_matrix_size + i * commit_rec->row_dim_bytes + CODEWORD_BYTES, 0);
  //     }
  //   }

  //   //Transpose block
  //   transpose_128_320(matrices_tmp.get() + commit_rec->transpose_matrix_size, matrices_tmp.get(), commit_rec->col_blocks);

  //   //The values matrix can be transposed using 8*128x128 matrix transpose. This matrix holds has values in the AES-1 position. The rest are all 0.
  //   matrix_buffer_values.EklundhBitTranspose(commit_rec->col_dim, commit_rec->row_dim_values);

  //   //Compute on block. Processes the block matrices in the same way as for the consistency check with the modification that there is no blinding values
  //   for (int l = 0; l < commit_rec->col_blocks; ++l) {
  //     for (int i = 0; i < CODEWORD_BITS; ++i) {
  //       //Load current row into val. If we are in one of the last AES_BITS commitments we directly add this to the final result as blinding. Else we multiply by alpha^(i+1) and store it in res_tmp.
  //       val = _mm_lddqu_si128((__m128i*) (matrices_tmp.get() + l * AES_BYTES + i * commit_rec->col_dim_bytes));

  //       if (j * commit_rec->col_dim + l * AES_BITS < num_values - AES_BITS) {
  //         //The actual commitments are multiplied with alpha
  //         mul128_karatsuba(val, alpha, &val_result[0], &val_result[1]);

  //         //Accumulate the val_result into res_tmp
  //         res_tmp[0][i] = _mm_xor_si128(res_tmp[0][i], val_result[0]);
  //         res_tmp[1][i] = _mm_xor_si128(res_tmp[1][i], val_result[1]);

  //         //Also, only process the values for the first AES_BITS bits as these are only this long
  //         if (i < AES_BITS) {
  //           val = _mm_lddqu_si128((__m128i*) (matrix_buffer_values.GetArr() + l * AES_BYTES + i * commit_rec->col_dim_bytes));
  //           mul128_karatsuba(val, alpha, &val_result[0], &val_result[1]);
  //           res_values_tmp[0][i] = _mm_xor_si128(res_values_tmp[0][i], val_result[0]);
  //           res_values_tmp[1][i] = _mm_xor_si128(res_values_tmp[1][i], val_result[1]);
  //         }
  //       } else if (j * commit_rec->col_dim + l * AES_BITS < num_values) {
  //         //The AES_BITS blinding one-time commitments are added directly to res_totals
  //         res_total[i] = _mm_xor_si128(res_total[i], val);
  //       }
  //     }
  //     //When done with one col_block we square the challenge element alpha. There are 8 col_blocks within each block
  //     gfmul128_no_refl(alpha, alpha, &alpha);
  //   }
  // }

  // //Need to manually delete this matrix
  // matrix_buffer_values.delCBitVector();

  // //mask is used to select the first SSEC linear combinations from res_total and store in final_result. Needed as we actually produce AES_BITS linear combinations due to convenience. However we only send and verify 2*SSEC of these.
  // uint8_t mask[CSEC_BYTES] = {0};
  // std::fill(mask, mask + SSEC_BYTES, 0xFF);
  // __m128i store_mask = _mm_lddqu_si128((__m128i*) mask);

  // //Reduce and store the resulting linear combinations
  // for (int i = 0; i < CODEWORD_BITS; ++i) {
  //   gfred128_no_refl(res_tmp[0][i], res_tmp[1][i], &res_tmp[0][i]);
  //   res_total[i] = _mm_xor_si128(res_total[i], res_tmp[0][i]);

  //   //Finally move the SSEC first linear combinations of the shares into final_result
  //   _mm_maskmoveu_si128(res_total[i], store_mask, (char*) (final_result + i * SSEC_BYTES));

  //   if (i < AES_BITS) {
  //     gfred128_no_refl(res_values_tmp[0][i], res_values_tmp[1][i], &res_values_tmp[0][i]);

  //     //Finally move the SSEC first linear combinations of the values into final_values_result
  //     _mm_maskmoveu_si128(res_values_tmp[0][i], store_mask, (char*) (final_values_result + i * SSEC_BYTES));
  //   }
  // }

  // //Receive the decommitments from CommitSnd that will be compared to the computed shares in final_result.
  // uint8_t decommit_shares0[2 * SSEC * CODEWORD_BYTES];
  // uint8_t* decommit_shares1 = decommit_shares0 + SSEC * CODEWORD_BYTES;
  // commit_rec->params.chan.ReceiveBlocking(decommit_shares0, 2 * SSEC * CODEWORD_BYTES);

  // //The below code is copied directly from the member function VerifyTransposedDecommits of CommitRec. Not the prettiest, but we need the values of matrix1 in order to perform our final consistency check
  // std::unique_ptr<uint8_t[]> matrix0(std::make_unique<uint8_t[]>(2 * commit_rec->transpose_matrix_size));
  // uint8_t* matrix1 = matrix0.get() + commit_rec->transpose_matrix_size;

  // //Read all shares in row-major order and compare choice bits row-wise. We check all num_values in parallel this way. Notice we select only the first num_values bits of each share using store_mask and bit-wise AND.
  // __m128i share, share0, share1, tmp;
  // for (int i = 0; i < CODEWORD_BITS; ++i) {
  //   share = _mm_lddqu_si128((__m128i*) (final_result + i * SSEC_BYTES));
  //   share = _mm_and_si128(share, store_mask);

  //   share0 = _mm_lddqu_si128((__m128i*) (decommit_shares0 + i * SSEC_BYTES));
  //   share0 = _mm_and_si128(share0, store_mask);

  //   share1 = _mm_lddqu_si128((__m128i*) (decommit_shares1 + i * SSEC_BYTES));
  //   share1 = _mm_and_si128(share1, store_mask);

  //   tmp = _mm_xor_si128(share0, share1);
  //   _mm_storeu_si128((__m128i*) (matrix0.get() + i * commit_rec->col_dim_single_bytes), tmp);

  //   if (GetBit(i, commit_rec->choices)) {
  //     if (!compare128(share, share1)) {
  //       return false;
  //     }
  //   } else {
  //     if (!compare128(share, share0)) {
  //       return false;
  //     }
  //   }
  // }
  // transpose_320_128(matrix0.get(), matrix1);

  // uint8_t tmp_checkbits[BCH_BYTES];
  // for (int i = 0; i < SSEC; ++i) {

  //   std::fill(tmp_checkbits, tmp_checkbits + BCH_BYTES, 0);

  //   commit_rec->code->Encode(matrix1 + i * commit_rec->row_dim_bytes, tmp_checkbits);

  //   if (!std::equal(tmp_checkbits, tmp_checkbits + BCH_BYTES, matrix1 + CSEC_BYTES + i * commit_rec->row_dim_bytes)) {
  //     std::cout << "Abort! Linear combination " << i << " is not a codeword" << std::endl;
  //     return false; //Not a codeword!
  //   }
  // }

  // //If the decommits check out, processed and verify that these were actually linear combinations of the postulated values as well. If everything is ok we are sure that the num_values postulated bits in values[] are indeed the committed lsbs and we return true.
  // for (int i = 0; i < SSEC; ++i) {
  //   if ((GetBitReversed(i, final_values_result + (AES_BYTES - 1) * SSEC)) != //the linear combinations of the values
  //       (GetLSB(matrix1 + i * commit_rec->row_dim_bytes) ^ //linear combinations of the decommitted values
  //        GetBit(i, delta_chal) ^ //Efficiently flips the bit if delta is included. Ensures that lsb(delta) = 1.
  //        GetBit(num_values - AES_BITS + i, values))) { //Masks out the blinding bit that is added for each combination
  //     return false;
  //   }
  // }

  return true;
}