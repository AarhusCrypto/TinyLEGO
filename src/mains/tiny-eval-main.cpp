#include "mains/mains.h"
#include "tiny/tiny-evaluator.h"

int main(int argc, const char* argv[]) {
  ezOptionParser opt;

  opt.overview = "TinyEvaluator Passing Parameters Guide.";
  opt.syntax = "Tinyeval first second third forth fifth sixth seventh";
  opt.example = "Tinyeval -n 4 -c aes -e 8,2,1 -o 0 -ip 10.11.100.216 -p 28001 -t 0\n\n";
  opt.footer = "ezOptionParser 0.1.4  Copyright (C) 2011 Remik Ziemlinski\nThis program is free and without warranty.\n";

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Display usage instructions.", // Help description.
    "-h",     // Flag token.
    "-help",  // Flag token.
    "--help", // Flag token.
    "--usage" // Flag token.
  );

  opt.add(
    default_num_iters.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Number of circuits to produce and evaluate.", // Help description.
    "-n"
  );

  opt.add(
    default_circuit_name.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Circuit name. Can be either aes, sha-1, sha-256 or cbc.", // Help description.
    "-c" // Flag token.
  );

  opt.add(
    default_execs.c_str(), // Default.
    0, // Required?
    3, // Number of args expected.
    ',', // Delimiter if expecting multiple args.
    "Number of parallel executions for each phase. Preprocessing, Offline and Online.", // Help description.
    "-e"
  );

  opt.add(
    default_optimize_online.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Optimize for online or overall efficiency", // Help description.
    "-o"
  );

  opt.add(
    default_ip_address.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "IP Address of Machine running TinyConst", // Help description.
    "-ip"
  );

  opt.add(
    default_port.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Port to listen on/connect to", // Help description.
    "-p"
  );

  opt.add(
    default_print_format.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Output timing info in special format", // Help description.
    "-t"
  );

  opt.add(
    default_circuit_file.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "circuit file", // Help description.
    "-f"
  );

  //Attempt to parse input
  opt.parse(argc, argv);

  //Check if help was requested and do some basic validation
  if (opt.isSet("-h")) {
    Usage(opt);
    return 1;
  }
  std::vector<std::string> badOptions;
  if (!opt.gotExpected(badOptions)) {
    for (int i = 0; i < badOptions.size(); ++i)
      std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i] << ".\n\n";
    Usage(opt);
    return 1;
  }

  //Copy inputs into the right variables
  int num_iters, pre_num_execs, offline_num_execs, online_num_execs, optimize_online, port, print_special_format;
  std::vector<int> num_execs;
  std::string circuit_name, circuit_file_name, ip_address;
  Circuit circuit;

  FILE* fileptr[2];
  uint8_t* buffer[2];
  uint8_t* expected_output;
  long filelen[2];
  opt.get("-n")->getInt(num_iters);
  opt.get("-c")->getString(circuit_name);
  opt.get("-f")->getString(circuit_file_name);

  opt.get("-e")->getInts(num_execs);
  pre_num_execs = num_execs[0];
  offline_num_execs = num_execs[1];
  online_num_execs = num_execs[2];

  opt.get("-o")->getInt(optimize_online);
  opt.get("-ip")->getString(ip_address);
  opt.get("-p")->getInt(port);
  opt.get("-t")->getInt(print_special_format);

  //Read the file_name_given
  if (!circuit_name.empty()) {
    if (!circuit_file_name.empty()) {
      std::cout << "Cannot both have -f and -c parameter set!" << std::endl;
      return 1;
    } else if (circuit_name.find("aes") != std::string::npos) {
      circuit_file_name = "test/data/AES-non-expanded.txt";
      fileptr[0] = fopen("test/data/aes_input_0.bin", "rb");
      fileptr[1] = fopen("test/data/aes_expected_0.bin", "rb");
    } else if (circuit_name.find("sha-256") != std::string::npos) {
      circuit_file_name = "test/data/sha-256.txt";
      fileptr[0] = fopen("test/data/sha256_input_0.bin", "rb");
      fileptr[1] = fopen("test/data/sha256_expected_0.bin", "rb");
    } else if (circuit_name.find("sha-1") != std::string::npos) {
      circuit_file_name = "test/data/sha-1.txt";
      fileptr[0] = fopen("test/data/sha1_input_0.bin", "rb");
      fileptr[1] = fopen("test/data/sha1_expected_0.bin", "rb");
    } else if (circuit_name.find("cbc") != std::string::npos) {
      circuit_file_name = "test/data/aescbcmac16.txt";
      fileptr[0] = fopen("test/data/cbc_input_0.bin", "rb");
      fileptr[1] = fopen("test/data/aes_expected_0.bin", "rb");
    } else {
      std::cout << "No circuit matching: " << circuit_name << ". Terminating" << std::endl;

      return 1;
    }

    //Read input and expected output from file to buffer[0] and buffer[1] and then set it to eval_input and expected output
    fseek(fileptr[0], 0, SEEK_END);
    fseek(fileptr[1], 0, SEEK_END);
    filelen[0] = ftell(fileptr[0]);
    filelen[1] = ftell(fileptr[1]);
    rewind(fileptr[0]);
    rewind(fileptr[1]);

    buffer[0] = new uint8_t[(filelen[0] + 1)];
    buffer[1] = new uint8_t[(filelen[1] + 1)];
    fread(buffer[0], filelen[0], 1, fileptr[0]);
    fread(buffer[1], filelen[1], 1, fileptr[1]);
  }

  if (circuit_file_name.empty()) {
    std::cout << "No circuit given" << std::endl;
    return 1;
  }

  circuit = read_text_circuit(circuit_file_name.c_str());

  std::unique_ptr<uint8_t[]> eval_input(std::make_unique<uint8_t[]>(BITS_TO_BYTES(circuit.num_eval_inp_wires))); //run on dummy 0 input

  //if predetermined case read actual input
  if (!circuit_name.empty()) {
    for (int i = 0; i < circuit.num_eval_inp_wires; ++i) {
      if (GetBitReversed(i, buffer[0]  + BITS_TO_BYTES(circuit.num_const_inp_wires))) {
        SetBit(i, 1, eval_input.get());
      }
    }
    expected_output = buffer[1];
  }

  //Compute number of gates, inputs and outputs that are to be preprocessed
  uint64_t num_gates = num_iters * circuit.num_and_gates;
  uint64_t num_inputs = num_iters * circuit.num_inp_wires;
  uint64_t num_outputs = num_iters * circuit.num_out_wires;

  std::vector<Circuit*> circuits;
  std::vector<uint8_t*> eval_inputs;
  std::vector<std::unique_ptr<uint8_t[]>> outputs;

  for (int i = 0; i < num_iters; ++i) {
    circuits.emplace_back(&circuit);
    eval_inputs.emplace_back(eval_input.get());
    outputs.emplace_back(std::make_unique<uint8_t[]>(BITS_TO_BYTES(circuit.num_out_wires)));
  }

  //Compute the required number of params that are to be created. We create one main param and one for each sub-thread that will be spawned later on. Need to know this at this point to setup context properly
  int num_params = std::max(pre_num_execs, offline_num_execs);
  num_params = std::max(num_params, online_num_execs);
  zmq::context_t context(NUM_IO_THREADS, 2 * (num_params + 1)); //We need two sockets pr. channel

  //Setup the main params object
  Params params(constant_seeds[1], num_gates, num_inputs, num_outputs, ip_address, (uint16_t) port, 1, context, pre_num_execs, GLOBAL_PARAMS_CHAN, optimize_online);

  TinyEvaluator tiny_eval(params);

  //Warm up network!
  uint8_t* dummy_val = new uint8_t[network_dummy_size]; //50 MB
  uint8_t* dummy_val_rec = new uint8_t[network_dummy_size]; //50 MB
  channel chan(OT_ADMIN_CHANNEL - 1, tiny_eval.ot_rec.net.rcvthread, tiny_eval.ot_rec.net.sndthread);
  chan.send(dummy_val, network_dummy_size);
  uint8_t* dummy_val_rec2 = chan.blocking_receive();

  params.chan.ReceiveBlocking(dummy_val_rec2, network_dummy_size);
  params.chan.Send(dummy_val, network_dummy_size);
  params.chan.bytes_received_vec[params.chan.received_pointer] = 0;
  params.chan.bytes_sent_vec[params.chan.sent_pointer] = 0;
  tiny_eval.ot_rec.net.m_vSocket->m_nSndCount = 0;
  tiny_eval.ot_rec.net.m_vSocket->m_nRcvCount = 0;
  delete[] dummy_val;
  delete[] dummy_val_rec;
  free(dummy_val_rec2);
  //Warm up network!

  //Run initial Setup (BaseOT) phase
  auto setup_begin = GET_TIME();
  mr_init_threading(); //Needed for Miracl library to work with threading.
  tiny_eval.Setup();
  mr_end_threading();
  auto setup_end = GET_TIME();

  uint64_t setup_sent = params.chan.GetTotalBytesSent() + tiny_eval.ot_rec.net.m_vSocket->getSndCnt();
  for (std::unique_ptr<Params>& thread_params : tiny_eval.thread_params_vec) {
    setup_sent += thread_params->chan.GetTotalBytesSent();
  }

  //Run Preprocessing phase
  auto preprocess_begin = GET_TIME();
  tiny_eval.Preprocess();
  auto preprocess_end = GET_TIME();

  uint64_t preprocess_sent = params.chan.GetTotalBytesSent() + tiny_eval.ot_rec.net.m_vSocket->getSndCnt();
  for (std::unique_ptr<Params>& thread_params : tiny_eval.thread_params_vec) {
    preprocess_sent += thread_params->chan.GetTotalBytesSent();
  }
  preprocess_sent -= setup_sent;

  //Preprocessing creates pre_num_execs sub-param objects. If more are needed in the offline and online phases we create them here.
  int extra_execs = num_params - params.num_execs;
  if (extra_execs < 0) {
  } else {
    std::unique_ptr<uint8_t[]> extra_thread_seeds(std::make_unique<uint8_t[]>(extra_execs * CSEC_BYTES));
    tiny_eval.params.rnd.GenRnd(extra_thread_seeds.get(), extra_execs * CSEC_BYTES);
    for (int i = 0; i < extra_execs; ++i) {
      tiny_eval.thread_params_vec.emplace_back(std::make_unique<Params>(tiny_eval.params, extra_thread_seeds.get() + i * CSEC_BYTES, tiny_eval.thread_params_vec[0]->num_pre_gates, tiny_eval.thread_params_vec[0]->num_pre_inputs, tiny_eval.thread_params_vec[0]->num_pre_outputs, params.num_execs + i));
      tiny_eval.commit_recs.emplace_back(std::make_unique<CommitReceiver>(*tiny_eval.thread_params_vec[params.num_execs + i], tiny_eval.rot_seeds.get(), tiny_eval.rot_choices.get()));
    }
  }

  // Figure out how many executions to run in offline phase
  int top_num_execs = std::min((int)circuits.size(), offline_num_execs);
  if (top_num_execs == 1) {
    tiny_eval.thread_pool.resize(top_num_execs);
  }
  // tiny_eval.thread_pool.resize(params.num_cpus * TP_MUL_FACTOR); //Very high performance benefit if on a high latency network as more executions can run in parallel!

  //Run Offline phase
  auto offline_begin = GET_TIME();
  tiny_eval.Offline(circuits, top_num_execs);
  auto offline_end = GET_TIME();

  uint64_t offline_sent = params.chan.GetTotalBytesSent() + tiny_eval.ot_rec.net.m_vSocket->getSndCnt();
  for (std::unique_ptr<Params>& thread_params : tiny_eval.thread_params_vec) {
    offline_sent += thread_params->chan.GetTotalBytesSent();
  }
  offline_sent -= (setup_sent + preprocess_sent);


  // If we are doing single evaluation then we have slightly better performance with a single thread running in the thread pool.
  int eval_num_execs = std::min((int)circuits.size(), online_num_execs);
  if (eval_num_execs == 1) {
    tiny_eval.thread_pool.resize(eval_num_execs);
  }

  //Prepare results vector. Needed as unique_ptr cannot be copied into Online call.
  std::vector<uint8_t*> outputs_raw;
  for (std::unique_ptr<uint8_t[]>& output : outputs) {
    outputs_raw.emplace_back(output.get());
  }

  //Run Online phase
  auto online_begin = GET_TIME();
  tiny_eval.Online(circuits, eval_inputs, outputs_raw, eval_num_execs);
  auto online_end = GET_TIME();

  uint64_t online_sent = params.chan.GetTotalBytesSent() + tiny_eval.ot_rec.net.m_vSocket->getSndCnt();
  for (std::unique_ptr<Params>& thread_params : tiny_eval.thread_params_vec) {
    online_sent += thread_params->chan.GetTotalBytesSent();
  }
  online_sent -= (setup_sent + preprocess_sent + offline_sent);

  //Check for correctness if predetermined case
  if (!circuit_name.empty()) {
    bool all_success = true;
    for (int i = 0; i < circuits.size(); ++i) {
      for (int j = 0; j < circuits[i]->num_out_wires; ++j) {
        if (GetBitReversed(j, expected_output) != GetBit(j, outputs[i].get())) {
          all_success = false;
        }
      }
    }
    if (!all_success && (circuit_file_name.find("test/data/aescbcmac16.txt") != std::string::npos)) { //Do not have expected results for AES-CBC-MAC circuit so we skip the check
      std::cout << "Wrong result!" << std::endl;
      //Do not return as we still report timings, even though result is wrong.
    }
  }

  // Average out the timings of each phase and print results
  uint64_t setup_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(setup_end - setup_begin).count();
  uint64_t preprocess_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(preprocess_end - preprocess_begin).count();
  uint64_t offline_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(offline_end - offline_begin).count();
  uint64_t online_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(online_end - online_begin).count();

  if (!print_special_format) {
    std::cout << "===== Eval timings for " << num_iters << " x " << circuit_file_name << "(" << (num_iters * circuit.num_and_gates) << ") with " << pre_num_execs << " preprocessing execs, " << top_num_execs << " offline execs and " << eval_num_execs << " online execs =====" << std::endl;

    std::cout << "Setup ms: " << (double) setup_time_nano / num_iters / 1000000 << ", data sent: " << (double) setup_sent / 1000 << " (" << (double) setup_sent / num_iters / 1000 << ")" << " kB" << std::endl;
    std::cout << "Preprocess ms: " << (double) preprocess_time_nano / num_iters / 1000000 << ", data sent: " << (double) preprocess_sent / 1000 << " (" << (double) preprocess_sent / num_iters / 1000 << ")" << " kB" << std::endl;
    std::cout << "Offline ms: " << (double) offline_time_nano / num_iters / 1000000 << ", data sent: " << (double) offline_sent / 1000 << " (" << (double) offline_sent / num_iters / 1000 << ")" << " kB" << std::endl;
    std::cout << "Online ms: " << (double) online_time_nano / num_iters / 1000000 << ", data sent: " << (double) online_sent / 1000 << " (" << (double) online_sent / num_iters / 1000 << ")" << " kB" << std::endl;
  } else {
    //Used for formatting output for paper
    std::cout << (double) setup_time_nano / num_iters / 1000000 << " " << (double) preprocess_time_nano / num_iters / 1000000 << " " << (double) offline_time_nano / num_iters / 1000000 << " " << (double) online_time_nano / num_iters / 1000000 << std::endl;
  }

  return 0;
}
