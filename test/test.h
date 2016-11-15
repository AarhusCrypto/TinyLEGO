#include "googletest/include/gtest/gtest.h"

//Hardcoded values for testing
static uint16_t default_port = 28001;
static std::string default_ip_address("localhost");
static int num_iters = 2;
static int test_num_gates = num_iters * 7000;
static int test_num_inputs = num_iters * 256;
static int test_num_outputs = num_iters * 128;