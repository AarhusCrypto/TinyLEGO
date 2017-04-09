#include "tiny-mains/mains.h"
#include "tiny/tiny-constructor.h"

int main(int argc, const char* argv[]) {
  ezOptionParser opt;

  opt.overview = "TinyConstructor Passing Parameters Guide.";
  opt.syntax = "Tinyconst first second third forth fifth sixth";
  opt.example = "Tinyconst -n 4 -c aes -e 8,2,1 -o 0 -ip 10.11.100.216 -p 28001\n\n";
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
  int num_iters, pre_num_execs, offline_num_execs, online_num_execs, optimize_online, port;
  std::vector<int> num_execs;
  std::string circuit_name, circuit_file_name, ip_address;
  Circuit circuit;

  FILE* fileptr;
  uint8_t* input_buffer;
  long filelen;
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

  //Read the file_name_given
  if (!circuit_name.empty()) {
    if (!circuit_file_name.empty()) {
      std::cout << "Cannot both have -f and -c parameter set!" << std::endl;
      return 1;
    } else if (circuit_name.find("aes") != std::string::npos) {
      circuit_file_name = "test/data/AES-non-expanded.txt";
      fileptr = fopen("test/data/aes_input_0.bin", "rb");
    } else if (circuit_name.find("sha-256") != std::string::npos) {
      circuit_file_name = "test/data/sha-256.txt";
      fileptr = fopen("test/data/sha256_input_0.bin", "rb");
    } else if (circuit_name.find("sha-1") != std::string::npos) {
      circuit_file_name = "test/data/sha-1.txt";
      fileptr = fopen("test/data/sha1_input_0.bin", "rb");
    } else if (circuit_name.find("cbc") != std::string::npos) {
      circuit_file_name = "test/data/aescbcmac16.txt";
      fileptr = fopen("test/data/cbc_input_0.bin", "rb");
    } else {

      std::cout << "No circuit matching: " << circuit_name << ". Terminating" << std::endl;
      return 1;
    }

    //Read input from file to input_buffer and then set it to const_input
    fseek(fileptr, 0, SEEK_END);
    filelen = ftell(fileptr);
    rewind(fileptr);
    input_buffer = new uint8_t[(filelen + 1)];
    fread(input_buffer, filelen, 1, fileptr);
  }

  if (circuit_file_name.empty()) {
    std::cout << "No circuit given" << std::endl;
    return 1;
  }

  circuit = read_text_circuit(circuit_file_name.c_str());

  osuCrypto::BitVector const_input(circuit.num_const_inp_wires);
  //if predetermined case read actual input
  if (!circuit_name.empty()) {
    for (int i = 0; i < circuit.num_const_inp_wires; ++i) {
      const_input[i] = GetBitReversed(i, input_buffer);
    }
  }

  //Compute number of gates, inputs and outputs that are to be preprocessed
  uint64_t num_gates = num_iters * circuit.num_and_gates;
  uint64_t num_inputs = num_iters * circuit.num_inp_wires;
  uint64_t num_outputs = num_iters * circuit.num_out_wires;

  std::vector<Circuit*> circuits;
  std::vector<osuCrypto::BitVector> const_inputs;
  for (int i = 0; i < num_iters; ++i) {
    circuits.emplace_back(&circuit);
    const_inputs.emplace_back(const_input);
  }

  //Compute the required number of params that are to be created. We create one main param and one for each sub-thread that will be spawned later on. Need to know this at this point to setup context properly
  int num_params = std::max(pre_num_execs, offline_num_execs);
  num_params = std::max(num_params, online_num_execs);

  //Setup the main params object
  Params params(num_gates, num_inputs, num_outputs, num_params);

  TinyConstructor tiny_const(tiny_constant_seeds[0], params);

  tiny_const.Connect(ip_address, (uint16_t) port);

  std::cout << "====== " << num_iters << " x " << circuit_name << " ======" << std::endl;

  //Values used for network syncing after each phase
  uint8_t rcv;
  uint8_t snd;

  //Run initial Setup (BaseOT) phase
  auto setup_begin = GET_TIME();
  tiny_const.Setup();
  auto setup_end = GET_TIME();

  uint64_t setup_data_sent = tiny_const.GetTotalDataSent();

  //Run Preprocessing phase
  auto preprocess_begin = GET_TIME();
  tiny_const.Preprocess();
  auto preprocess_end = GET_TIME();

  //Sync with Evaluator
  tiny_const.chan.recv(&rcv, 1);
  tiny_const.chan.send(&snd, 1);

  uint64_t preprocess_data_sent = tiny_const.GetTotalDataSent() - setup_data_sent;

  // Figure out how many executions to run in offline phase
  int top_num_execs = std::min((int)circuits.size(), offline_num_execs);
  if (top_num_execs == 1) {
    tiny_const.thread_pool.resize(top_num_execs);
  }

  // tiny_const.thread_pool.resize(params.num_cpus * TP_MUL_FACTOR); //Very high performance benefit if on a high latency network as more executions can run in parallel!

  //Run Offline phase
  auto offline_begin = GET_TIME();
  tiny_const.Offline(circuits, top_num_execs);
  auto offline_end = GET_TIME();

  //Sync with Evaluator
  tiny_const.chan.recv(&rcv, 1);
  tiny_const.chan.send(&snd, 1);

  uint64_t offline_data_sent = tiny_const.GetTotalDataSent() - setup_data_sent - preprocess_data_sent;

  // If we are doing single evaluation then we have slightly better performance with a single thread running in the thread pool.
  int eval_num_execs = std::min((int)circuits.size(), online_num_execs);
  if (eval_num_execs == 1) {
    tiny_const.thread_pool.resize(eval_num_execs);
  }

  //Run Online phase
  auto online_begin = GET_TIME();
  tiny_const.Online(circuits, const_inputs, eval_num_execs);
  auto online_end = GET_TIME();

  //Sync with Evaluator
  tiny_const.chan.recv(&rcv, 1);
  tiny_const.chan.send(&snd, 1);

  uint64_t online_data_sent = tiny_const.GetTotalDataSent() - setup_data_sent - preprocess_data_sent - offline_data_sent;

  // Average out the timings of each phase and print results
  uint64_t setup_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(setup_end - setup_begin).count();
  uint64_t preprocess_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(preprocess_end - preprocess_begin).count();
  uint64_t offline_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(offline_end - offline_begin).count();
  uint64_t online_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(online_end - online_begin).count();

  std::cout << "Setup ms: " << (double) setup_time_nano / num_iters / 1000000 << ", data sent: " << (double) setup_data_sent / num_iters / 1000 << " kB" << std::endl;
  std::cout << "Preprocess ms: " << (double) preprocess_time_nano / num_iters / 1000000 << ", data sent: " << (double) preprocess_data_sent / num_iters / 1000 << " kB" << std::endl;
  std::cout << "Offline ms: " << (double) offline_time_nano / num_iters / 1000000 << ", data sent: " << (double) offline_data_sent / num_iters / 1000 << " kB" << std::endl;
  std::cout << "Online ms: " << (double) online_time_nano / num_iters / 1000000 << ", data sent: " << (double) online_data_sent / num_iters / 1000 << " kB" << std::endl;

  return 0;
}