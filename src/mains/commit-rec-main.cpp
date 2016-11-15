#include "mains/mains.h"
#include "tiny/tiny-evaluator.h"

int main(int argc, const char* argv[]) {
  ezOptionParser opt;

  opt.overview = "Commit Receiver Passing Parameters Guide.";
  opt.syntax = "Commitrec first second";
  opt.example = "Commitrec -n 10000 -e 4 -ip 10.11.100.216 -p 28001 -t 0\n\n";
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
    default_num_commits.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Number of commits to produce and decommit.", // Help description.
    "-n"
  );

  opt.add(
    default_num_commit_execs.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Number of parallel executions to run. These will share the workload.", // Help description.
    "-e" // Flag token.
  );

  opt.add(
    default_ip_address.c_str(), // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "IP Address of Machine running Commitsnd", // Help description.
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
  int num_commits, num_execs, port, print_special_format;
  std::string ip_address;

  opt.get("-n")->getInt(num_commits);
  opt.get("-e")->getInt(num_execs);
  opt.get("-p")->getInt(port);
  opt.get("-ip")->getString(ip_address);
  opt.get("-t")->getInt(print_special_format);

  zmq::context_t context(NUM_IO_THREADS, 2 * (num_execs + 1));
  Params params(constant_seeds[1], commit_dummy_size, commit_dummy_input, commit_dummy_output, ip_address, (uint16_t) port, 1, context, num_execs); //Hardcoded dummy values
  ctpl::thread_pool thread_pool(params.num_cpus);
  params.num_commits = num_commits;
  params.num_OT = CODEWORD_BITS;

  ALSZDOTExtRec ot_rec(params);
  

  //Warm up network!
  uint8_t* dummy_val = new uint8_t[network_dummy_size]; //50 MB
  uint8_t* dummy_val_rec = new uint8_t[network_dummy_size]; //50 MB
  channel chan(OT_ADMIN_CHANNEL - 1, ot_rec.net.rcvthread, ot_rec.net.sndthread);
  chan.send(dummy_val, network_dummy_size);
  uint8_t* dummy_val_rec2 = chan.blocking_receive();
  params.chan.ReceiveBlocking(dummy_val_rec2, network_dummy_size);
  params.chan.Send(dummy_val, network_dummy_size);
  params.chan.bytes_received_vec[params.chan.received_pointer] = 0;
  params.chan.bytes_sent_vec[params.chan.sent_pointer] = 0;
  delete[] dummy_val;
  delete[] dummy_val_rec;
  free(dummy_val_rec2);
  //Warm up network!

  //Base OTs
  auto baseot_begin = GET_TIME();
  mr_init_threading();
  ot_rec.InitOTReceiver();
  mr_end_threading();
  auto baseot_end = GET_TIME();

  auto otx_begin = GET_TIME();
  ot_rec.Receive();
  auto otx_end = GET_TIME();

  std::unique_ptr<uint8_t[]> rot_seeds(std::make_unique<uint8_t[]>(CODEWORD_BITS * CSEC_BYTES));
  std::unique_ptr<uint8_t[]> rot_choices(std::make_unique<uint8_t[]>(BITS_TO_BYTES(CODEWORD_BITS)));

  int rot_start_pos = params.num_OT - CODEWORD_BITS;
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    if (GetBit(rot_start_pos + i, ot_rec.choices_outer.get())) {
      SetBit(i, 1, rot_choices.get());
    } else {
      SetBit(i, 0, rot_choices.get());
    }
  }

  //Hash away the global correlation and store in rot_seeds
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    params.crypt.hash(rot_seeds.get() + i * CSEC_BYTES, CSEC_BYTES, ot_rec.response_outer.get() + (rot_start_pos + i) * CSEC_BYTES, CSEC_BYTES);
  }

  auto ot_end = GET_TIME();

  auto commit_begin = GET_TIME();

  std::vector<std::future<void>> execs_finished(params.num_execs);
  std::unique_ptr<uint8_t[]> thread_seeds(std::make_unique<uint8_t[]>(CSEC_BYTES * params.num_execs));
  params.rnd.GenRnd(thread_seeds.get(), CSEC_BYTES * params.num_execs);

  std::vector<std::unique_ptr<Params>> thread_params_vec;
  std::vector<std::unique_ptr<CommitReceiver>> commit_recs;

  for (int exec_id = 0; exec_id < params.num_execs; ++exec_id) {
    thread_params_vec.emplace_back(std::make_unique<Params>(thread_seeds.get() + exec_id * CSEC_BYTES, params.num_pre_gates, params.num_pre_inputs, params.num_pre_outputs, params.ip_address, params.port, params.net_role, params.context, params.num_execs, exec_id));
    Params* thread_params = thread_params_vec[exec_id].get();
    thread_params->num_commits = params.num_commits / params.num_execs;
    thread_params->num_OT = CODEWORD_BITS;

    commit_recs.emplace_back(std::make_unique<CommitReceiver>(*thread_params, rot_seeds.get(), rot_choices.get()));
    CommitReceiver* commit_rec = commit_recs[exec_id].get();

    execs_finished[exec_id] = thread_pool.push([thread_params, commit_rec, exec_id] (int id) {
      commit_rec->Commit();
    });
  }

  for (std::future<void>& r : execs_finished) {
    r.wait();
  }

  auto commit_end = GET_TIME();

  std::vector<std::future<void>> decommits_finished(params.num_execs);

  auto decommit_begin = GET_TIME();
  for (int exec_id = 0; exec_id < params.num_execs; ++exec_id) {
    Params* thread_params = thread_params_vec[exec_id].get();

    //We store our local state in the containers as we need to access them for future use
    CommitReceiver* commit_rec = commit_recs[exec_id].get();

    decommits_finished[exec_id] = thread_pool.push([thread_params, commit_rec, exec_id] (int id) {
      //Decommit
      std::unique_ptr<uint8_t[]> data(std::make_unique<uint8_t[]>(thread_params->num_commits * (CODEWORD_BYTES + CSEC_BYTES)));
      uint8_t* computed_shares = data.get();
      uint8_t* values = computed_shares + thread_params->num_commits * CODEWORD_BYTES;

      for (int i = 0; i < thread_params->num_commits; ++i) {
        std::copy(commit_rec->commit_shares[i], commit_rec->commit_shares[i] + CODEWORD_BYTES, computed_shares + i * CODEWORD_BYTES);
      }
      thread_params->chan.ReceiveBlocking(values, thread_params->num_commits * CSEC_BYTES);
      if (!commit_rec->BatchDecommit(computed_shares, thread_params->num_commits, values)) {
        std::cout << "Error in Decommit!" << std::endl;
      }
    });
  }

  for (std::future<void>& r : decommits_finished) {
    r.wait();
  }
  auto decommit_end = GET_TIME();


  uint64_t base_ot_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(baseot_end - baseot_begin).count();
  uint64_t otx_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(otx_end - otx_begin).count();
  uint64_t commit_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(commit_end - commit_begin).count();
  uint64_t decommit_time_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(decommit_end - decommit_begin).count();

  if (!print_special_format) {

    std::cout << "===== Timings for receiver doing " << num_commits << " random commits using " << num_execs << " parallel execs " << std::endl;

    std::cout << "OT ms: " << (double) base_ot_time_nano / 1000000 << std::endl;
    std::cout << "OTX ms: " << (double) otx_time_nano / 1000000 << std::endl;
    std::cout << "Amortized OT ms: " << (double) (base_ot_time_nano + otx_time_nano) / params.num_commits / 1000000 << std::endl;
    std::cout << "Commit us (with OT): " << (double) (commit_time_nano + base_ot_time_nano) / params.num_commits / 1000 << std::endl;
    std::cout << "Commit us: " << (double) commit_time_nano / params.num_commits / 1000 << std::endl;
    std::cout << "Commit total ms: " << (double) (commit_time_nano + base_ot_time_nano) / 1000000 << std::endl;
    std::cout << "Decommit us: " << (double) decommit_time_nano / params.num_commits / 1000 << std::endl;
    std::cout << "Decommit total ms: " << (double) decommit_time_nano / 1000000 << std::endl;
  } else {
    std::cout << params.num_commits << " " << num_execs << " " << (double) commit_time_nano / params.num_commits / 1000 << " " << (double) (commit_time_nano + base_ot_time_nano) / params.num_commits / 1000 << " " << (double) decommit_time_nano / params.num_commits / 1000 << std::endl;
  }

  return 0;
}
