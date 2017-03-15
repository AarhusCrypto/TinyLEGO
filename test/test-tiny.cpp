#include "test.h"

#include "tiny/tiny-constructor.h"
#include "tiny/tiny-evaluator.h"
#include "tiny/params.h"

void RunConst(TinyConstructor& tiny_const, std::vector<Circuit*>& circuits, std::vector<uint8_t*>& inputs) {

  tiny_const.Connect(default_ip_address, default_port);
  tiny_const.Setup();
  tiny_const.Preprocess();
  tiny_const.Offline(circuits, tiny_const.params.num_execs);
  tiny_const.Online(circuits, inputs, tiny_const.params.num_execs);

}

void RunEval(TinyEvaluator& tiny_eval, std::vector<Circuit*>& circuits, std::vector<uint8_t*>& inputs, std::vector<std::unique_ptr<uint8_t[]>>& outputs, uint8_t* expected_output) {

  std::vector<uint8_t*> outputs_raw;
  for(std::unique_ptr<uint8_t[]>& ptr: outputs) {
    outputs_raw.emplace_back(ptr.get());
  }

  tiny_eval.Connect(default_ip_address, default_port);
  tiny_eval.Setup();
  tiny_eval.Preprocess();
  tiny_eval.Offline(circuits, tiny_eval.params.num_execs);
  tiny_eval.Online(circuits, inputs, outputs_raw, tiny_eval.params.num_execs);

  // for (int i = 0; i < circuits.size(); ++i) {
  //   for (int j = 0; j < circuits[i]->num_out_wires; ++j) {

  //     ASSERT_TRUE(GetBitReversed(j, expected_output) == GetBit(j, outputs[i].get()));
  //   }
  // }
}

// TEST(Protocol, AES) {
//   zmq::context_t context0(1);
//   zmq::context_t context1(1);
//   Params params_const(constant_seeds[0], num_iters * 7000, num_iters * 256, num_iters * 128, default_ip_address, default_port, 0, context0, 2, GLOBAL_PARAMS_CHAN);
//   Params params_eval(constant_seeds[1],  num_iters * 7000, num_iters * 256, num_iters * 128, default_ip_address, default_port, 1, context1, 2, GLOBAL_PARAMS_CHAN);

//   TinyConstructor tiny_const(params_const);
//   TinyEvaluator tiny_eval(params_eval);
//   size_t size;
//   Circuit circuit = read_text_circuit("test/data/AES-non-expanded.txt");

//   FILE *fileptr[2];
//   fileptr[0] = fopen("test/data/aes_input_0.bin", "rb");  // Open the file in binary mode
//   fileptr[1] = fopen("test/data/aes_expected_0.bin", "rb");  // Open the file in binary mode
//   uint8_t *buffer[2];
//   long filelen[2];
//   fseek(fileptr[0], 0, SEEK_END);          // Jump to the end of the file
//   fseek(fileptr[1], 0, SEEK_END);          // Jump to the end of the file
//   filelen[0] = ftell(fileptr[0]);             // Get the current byte offset in the file
//   filelen[1] = ftell(fileptr[1]);             // Get the current byte offset in the file
//   rewind(fileptr[0]);                      // Jump back to the beginning of the file
//   rewind(fileptr[1]);                      // Jump back to the beginning of the file

//   buffer[0] = new uint8_t[(filelen[0] + 1)]; // Enough memory for file + \0
//   buffer[1] = new uint8_t[(filelen[1] + 1)]; // Enough memory for file + \0
//   fread(buffer[0], filelen[0], 1, fileptr[0]); // Read in the entire file
//   fread(buffer[1], filelen[1], 1, fileptr[1]); // Read in the entire file

//   uint8_t* const_input = new uint8_t[BITS_TO_BYTES(circuit.num_const_inp_wires)];
//   //Read input the right way!
//   for (int i = 0; i < circuit.num_const_inp_wires; ++i) {
//     if (GetBitReversed(i, buffer[0])) {
//       SetBit(i, 1, const_input);
//     } else {
//       SetBit(i, 0, const_input);
//     }
//   }

//   uint8_t* eval_input = new uint8_t[BITS_TO_BYTES(circuit.num_eval_inp_wires)];
//   for (int i = 0; i < circuit.num_eval_inp_wires; ++i) {
//     if (GetBitReversed(i, buffer[0] +  + BITS_TO_BYTES(circuit.num_const_inp_wires))) {
//       SetBit(i, 1, eval_input);
//     } else {
//       SetBit(i, 0, eval_input);
//     }
//   }

//   uint8_t* expected_output = buffer[1];

//   std::vector<Circuit*> circuits;
//   std::vector<uint8_t*> eval_inputs;
//   std::vector<uint8_t*> const_inputs;
//   std::vector<std::unique_ptr<uint8_t[]>> outputs;
//   for (int i = 0; i < num_iters; ++i) {
//     circuits.emplace_back(&circuit);
//     const_inputs.emplace_back(const_input);
//     eval_inputs.emplace_back(eval_input);
//     outputs.emplace_back(new uint8_t[BITS_TO_BYTES(circuit.num_out_wires)]);
//   }

//   mr_init_threading();
//   thread tiny_const_thread(RunConst, std::ref(tiny_const), std::ref(circuits), std::ref(const_inputs));
//   thread tiny_eval_thread(RunEval, std::ref(tiny_eval), std::ref(circuits), std::ref(eval_inputs), std::ref(outputs), expected_output);

//   tiny_const_thread.join();
//   tiny_eval_thread.join();
//   mr_end_threading();
// }

TEST(Protocol, AES) {
  Params params_const(num_iters * 7000, num_iters * 256, num_iters * 128, 2);
  Params params_eval(num_iters * 7000, num_iters * 256, num_iters * 128, 2);
  
  TinyConstructor tiny_const(tiny_constant_seeds[0], params_const);
  TinyEvaluator tiny_eval(tiny_constant_seeds[1], params_eval);
  
  size_t size;
  Circuit circuit = read_text_circuit("test/data/AES-non-expanded.txt");

  FILE *fileptr[2];
  fileptr[0] = fopen("test/data/aes_input_0.bin", "rb");  // Open the file in binary mode
  fileptr[1] = fopen("test/data/aes_expected_0.bin", "rb");  // Open the file in binary mode
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

  uint8_t* const_input = new uint8_t[BITS_TO_BYTES(circuit.num_const_inp_wires)];
  //Read input the right way!
  for (int i = 0; i < circuit.num_const_inp_wires; ++i) {
    if (GetBitReversed(i, buffer[0])) {
      SetBit(i, 1, const_input);
    } else {
      SetBit(i, 0, const_input);
    }
  }

  uint8_t* eval_input = new uint8_t[BITS_TO_BYTES(circuit.num_eval_inp_wires)];
  for (int i = 0; i < circuit.num_eval_inp_wires; ++i) {
    if (GetBitReversed(i, buffer[0] +  + BITS_TO_BYTES(circuit.num_const_inp_wires))) {
      SetBit(i, 1, eval_input);
    } else {
      SetBit(i, 0, eval_input);
    }
  }

  uint8_t* expected_output = buffer[1];

  std::vector<Circuit*> circuits;
  std::vector<uint8_t*> eval_inputs;
  std::vector<uint8_t*> const_inputs;
  std::vector<std::unique_ptr<uint8_t[]>> outputs;
  for (int i = 0; i < num_iters; ++i) {
    circuits.emplace_back(&circuit);
    const_inputs.emplace_back(const_input);
    eval_inputs.emplace_back(eval_input);
    outputs.emplace_back(new uint8_t[BITS_TO_BYTES(circuit.num_out_wires)]);
  }

  std::thread tiny_const_thread(RunConst, std::ref(tiny_const), std::ref(circuits), std::ref(const_inputs));
  std::thread tiny_eval_thread(RunEval, std::ref(tiny_eval), std::ref(circuits), std::ref(eval_inputs), std::ref(outputs), expected_output);
  
  tiny_const_thread.join();
  tiny_eval_thread.join();
}

// TEST(Protocol, DISABLED_SHA_1) {
//   num_iters = 2;
//   zmq::context_t context0(1);
//   zmq::context_t context1(1);
//   Params params_const(constant_seeds[0], num_iters * 236624, num_iters * 512, num_iters * 256, default_ip_address, default_port+2000, 0, context0, 1);
//   Params params_eval(constant_seeds[1],  num_iters * 236624, num_iters * 512, num_iters * 256, default_ip_address, default_port+2000, 1, context1, 1);

//   TinyConstructor tiny_const(params_const);
//   TinyEvaluator tiny_eval(params_eval);

//   size_t size;
//   Circuit circuit = read_text_circuit("test/data/sha-1.txt");

//   FILE *fileptr[2];
//   fileptr[0] = fopen("test/data/sha1_input_0.bin", "rb");  // Open the file in binary mode
//   fileptr[1] = fopen("test/data/sha1_expected_0.bin", "rb");  // Open the file in binary mode
//   uint8_t *buffer[2];
//   long filelen[2];
//   fseek(fileptr[0], 0, SEEK_END);          // Jump to the end of the file
//   fseek(fileptr[1], 0, SEEK_END);          // Jump to the end of the file
//   filelen[0] = ftell(fileptr[0]);             // Get the current byte offset in the file
//   filelen[1] = ftell(fileptr[1]);             // Get the current byte offset in the file
//   rewind(fileptr[0]);                      // Jump back to the beginning of the file
//   rewind(fileptr[1]);                      // Jump back to the beginning of the file

//   buffer[0] = new uint8_t[(filelen[0] + 1)]; // Enough memory for file + \0
//   buffer[1] = new uint8_t[(filelen[1] + 1)]; // Enough memory for file + \0
//   fread(buffer[0], filelen[0], 1, fileptr[0]); // Read in the entire file
//   fread(buffer[1], filelen[1], 1, fileptr[1]); // Read in the entire file

//   uint8_t* const_input = new uint8_t[BITS_TO_BYTES(circuit.num_const_inp_wires)];
//   //Read input the right way!
//   for (int i = 0; i < circuit.num_const_inp_wires; ++i) {
//     if (GetBitReversed(i, buffer[0])) {
//       SetBit(i, 1, const_input);
//     } else {
//       SetBit(i, 0, const_input);
//     }
//   }

//   uint8_t* eval_input = NULL;
//   uint8_t* expected_output = buffer[1];

//   std::vector<Circuit*> circuits;
//   std::vector<uint8_t*> eval_inputs;
//   std::vector<uint8_t*> const_inputs;
//   std::vector<std::unique_ptr<uint8_t[]>> outputs;
//   for (int i = 0; i < num_iters; ++i) {
//     circuits.emplace_back(&circuit);
//     const_inputs.emplace_back(const_input);
//     eval_inputs.emplace_back(eval_input);
//     outputs.emplace_back(new uint8_t[BITS_TO_BYTES(circuit.num_out_wires)]);
//   }


//   mr_init_threading();
//   thread tiny_const_thread(RunConst, std::ref(tiny_const), std::ref(circuits), std::ref(const_inputs));
//   thread tiny_eval_thread(RunEval, std::ref(tiny_eval), std::ref(circuits), std::ref(eval_inputs), std::ref(outputs), expected_output);

//   tiny_const_thread.join();
//   tiny_eval_thread.join();
//   mr_end_threading();
// }

// TEST(Protocol, DISABLED_SHA_256) {
//   num_iters = 2;
//   zmq::context_t context0(1);
//   zmq::context_t context1(1);
//   Params params_const(constant_seeds[0], num_iters * 236624, num_iters * 512, num_iters * 256, default_ip_address, default_port+2000, 0, context0, 1);
//   Params params_eval(constant_seeds[1],  num_iters * 236624, num_iters * 512, num_iters * 256, default_ip_address, default_port+2000, 1, context1, 1);

//   TinyConstructor tiny_const(params_const);
//   TinyEvaluator tiny_eval(params_eval);

//   size_t size;
//   Circuit circuit = read_text_circuit("test/data/sha-256.txt");

//   FILE *fileptr[2];
//   fileptr[0] = fopen("test/data/sha256_input_0.bin", "rb");  // Open the file in binary mode
//   fileptr[1] = fopen("test/data/sha256_expected_0.bin", "rb");  // Open the file in binary mode
//   uint8_t *buffer[2];
//   long filelen[2];
//   fseek(fileptr[0], 0, SEEK_END);          // Jump to the end of the file
//   fseek(fileptr[1], 0, SEEK_END);          // Jump to the end of the file
//   filelen[0] = ftell(fileptr[0]);             // Get the current byte offset in the file
//   filelen[1] = ftell(fileptr[1]);             // Get the current byte offset in the file
//   rewind(fileptr[0]);                      // Jump back to the beginning of the file
//   rewind(fileptr[1]);                      // Jump back to the beginning of the file

//   buffer[0] = new uint8_t[(filelen[0] + 1)]; // Enough memory for file + \0
//   buffer[1] = new uint8_t[(filelen[1] + 1)]; // Enough memory for file + \0
//   fread(buffer[0], filelen[0], 1, fileptr[0]); // Read in the entire file
//   fread(buffer[1], filelen[1], 1, fileptr[1]); // Read in the entire file

//   uint8_t* const_input = new uint8_t[BITS_TO_BYTES(circuit.num_const_inp_wires)];
//   //Read input the right way!
//   for (int i = 0; i < circuit.num_const_inp_wires; ++i) {
//     if (GetBitReversed(i, buffer[0])) {
//       SetBit(i, 1, const_input);
//     } else {
//       SetBit(i, 0, const_input);
//     }
//   }

//   uint8_t* eval_input = NULL;
//   uint8_t* expected_output = buffer[1];

//   std::vector<Circuit*> circuits;
//   std::vector<uint8_t*> eval_inputs;
//   std::vector<uint8_t*> const_inputs;
//   std::vector<std::unique_ptr<uint8_t[]>> outputs;
//   for (int i = 0; i < num_iters; ++i) {
//     circuits.emplace_back(&circuit);
//     const_inputs.emplace_back(const_input);
//     eval_inputs.emplace_back(eval_input);
//     outputs.emplace_back(new uint8_t[BITS_TO_BYTES(circuit.num_out_wires)]);
//   }


//   mr_init_threading();
//   thread tiny_const_thread(RunConst, std::ref(tiny_const), std::ref(circuits), std::ref(const_inputs));
//   thread tiny_eval_thread(RunEval, std::ref(tiny_eval), std::ref(circuits), std::ref(eval_inputs), std::ref(outputs), expected_output);

//   tiny_const_thread.join();
//   tiny_eval_thread.join();
//   mr_end_threading();
// }