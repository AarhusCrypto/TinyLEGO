#include "ezOptionParser/ezOptionParser.hpp"

using namespace ez;

//Hardcoded constants
static const uint32_t network_dummy_size = 50000000;
static const uint32_t commit_dummy_size = 1000000;
static const uint32_t commit_dummy_input = 128;
static const uint32_t commit_dummy_output = 128;

//Hardcoded default values
static std::string default_num_iters("10");
static std::string default_circuit_name("");
static std::string default_circuit_file("");
static std::string default_execs("1, 1, 1");
static std::string default_optimize_online("0");
static std::string default_ip_address("localhost");
static std::string default_port("28001");
static std::string default_print_format("0");

static std::string default_num_commits("10000");
static std::string default_num_commit_execs("1");

void Usage(ezOptionParser& opt) {
  std::string usage;
  opt.getUsage(usage);
  std::cout << usage;
};