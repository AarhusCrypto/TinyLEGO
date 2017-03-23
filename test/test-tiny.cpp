#include "test.h"

#include "tiny/tiny-constructor.h"
#include "tiny/tiny-evaluator.h"
#include "tiny/params.h"

void RunConst(TinyConstructor& tiny_const, std::vector<Circuit*>& circuits, std::vector<osuCrypto::BitVector>& inputs) {

  tiny_const.Connect(default_ip_address, default_port);
  tiny_const.Setup();
  tiny_const.Preprocess();
  tiny_const.Offline(circuits, tiny_const.params.num_max_execs);
  tiny_const.Online(circuits, inputs, tiny_const.params.num_max_execs);

}

void RunEval(TinyEvaluator& tiny_eval, std::vector<Circuit*>& circuits, std::vector<osuCrypto::BitVector>& inputs, osuCrypto::BitVector& expected_output) {

  std::vector<osuCrypto::BitVector> outputs_raw;
  for (int i = 0; i < circuits.size(); ++i) {
    outputs_raw.emplace_back(circuits[i]->num_out_wires);
  }

  tiny_eval.Connect(default_ip_address, default_port);
  tiny_eval.Setup();
  tiny_eval.Preprocess();
  tiny_eval.Offline(circuits, tiny_eval.params.num_max_execs);
  tiny_eval.Online(circuits, inputs, outputs_raw, tiny_eval.params.num_max_execs);

  for (int i = 0; i < circuits.size(); ++i) {
    for (int j = 0; j < circuits[i]->num_out_wires; ++j) {

      ASSERT_TRUE(GetBitReversed(j, expected_output.data()) == outputs_raw[i][j]);
    }
  }
}

void TestCircuit(std::string circuit_file, std::string input_file, std::string exp_output_file, uint32_t num_iters) {

  Circuit circuit = read_text_circuit(circuit_file.c_str());

  FILE *fileptr[2];
  fileptr[0] = fopen(input_file.c_str(), "rb");  // Open the file in binary mode
  fileptr[1] = fopen(exp_output_file.c_str(), "rb");  // Open the file in binary mode
  uint8_t *buffer[2];
  long filelen[2];
  fseek(fileptr[0], 0, SEEK_END);          // Jump to the end of the file
  fseek(fileptr[1], 0, SEEK_END);          // Jump to the end of the file
  filelen[0] = ftell(fileptr[0]);             // Get the current byte offset in the file
  filelen[1] = ftell(fileptr[1]);             // Get the current byte offset in the file
  rewind(fileptr[0]);                      // Jump back to the beginning of the file
  rewind(fileptr[1]);                      // Jump back to the beginning of the file

  buffer[0] = new uint8_t[(filelen[0] + 1)]; // Enough memory for file + \0
  buffer[1] = new uint8_t[(filelen[1] + 1)]; // Enough memory for file + \0
  fread(buffer[0], filelen[0], 1, fileptr[0]); // Read in the entire file
  fread(buffer[1], filelen[1], 1, fileptr[1]); // Read in the entire file

  Params params_const(num_iters * circuit.num_and_gates, num_iters * circuit.num_inp_wires, num_iters * circuit.num_out_wires, 2);
  Params params_eval(num_iters * circuit.num_and_gates, num_iters * circuit.num_inp_wires, num_iters * circuit.num_out_wires, 2);

  TinyConstructor tiny_const(tiny_constant_seeds[0], params_const);
  TinyEvaluator tiny_eval(tiny_constant_seeds[1], params_eval);

  //Read input the right way!

  osuCrypto::BitVector const_input(circuit.num_const_inp_wires);
  for (int i = 0; i < circuit.num_const_inp_wires; ++i) {
    const_input[i] = GetBitReversed(i, buffer[0]);
  }

  osuCrypto::BitVector eval_input(circuit.num_eval_inp_wires);
  for (int i = 0; i < circuit.num_eval_inp_wires; ++i) {
    eval_input[i] = GetBitReversed(i, buffer[0] + BITS_TO_BYTES(circuit.num_const_inp_wires));
  }

  osuCrypto::BitVector expected_outputs(buffer[1], circuit.num_out_wires);

  std::vector<Circuit*> circuits;
  std::vector<osuCrypto::BitVector> eval_inputs;
  std::vector<osuCrypto::BitVector> const_inputs;
  for (int i = 0; i < num_iters; ++i) {
    circuits.emplace_back(&circuit);
    const_inputs.emplace_back(const_input);
    eval_inputs.emplace_back(eval_input);
  }

  std::thread tiny_const_thread(RunConst, std::ref(tiny_const), std::ref(circuits), std::ref(const_inputs));
  std::thread tiny_eval_thread(RunEval, std::ref(tiny_eval), std::ref(circuits), std::ref(eval_inputs), std::ref(expected_outputs));

  tiny_const_thread.join();
  tiny_eval_thread.join();
}

TEST(Protocol, AES) {

  TestCircuit("test/data/AES-non-expanded.txt", "test/data/aes_input_0.bin", "test/data/aes_expected_0.bin", num_iters);
}

TEST(Protocol, SHA_1) {

  TestCircuit("test/data/sha-1.txt", "test/data/sha1_input_0.bin", "test/data/sha1_expected_0.bin", 2);
}

TEST(Protocol, SHA_256) {

  TestCircuit("test/data/sha-256.txt", "test/data/sha256_input_0.bin", "test/data/sha256_expected_0.bin", 2);
}