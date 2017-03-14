#include "test.h"

#include "circuit/circuit.h"
#include "tiny-util/util.h"

TEST(GetCircuit, Parse) {
  Circuit c = read_text_circuit("test/data/AES-non-expanded.txt");

  FILE *fileptr[2];
  uint8_t *buffer[2];
  long filelen[2];

  fileptr[0] = fopen("test/data/aes_input_0.bin", "rb");  // Open the file in binary mode
  fileptr[1] = fopen("test/data/aes_expected_0.bin", "rb");  // Open the file in binary mode
  fseek(fileptr[0], 0, SEEK_END);          // Jump to the end of the file
  fseek(fileptr[1], 0, SEEK_END);          // Jump to the end of the file
  filelen[0] = ftell(fileptr[0]);             // Get the current byte offset in the file
  filelen[1] = ftell(fileptr[1]);             // Get the current byte offset in the file
  rewind(fileptr[0]);                      // Jump back to the beginning of the file
  rewind(fileptr[1]);                      // Jump back to the beginning of the file

  buffer[0] = new uint8_t[(filelen[0] + 1) * sizeof(uint8_t)]; // Enough memory for file + \0
  buffer[1] = new uint8_t[(filelen[1] + 1) * sizeof(uint8_t)]; // Enough memory for file + \0
  fread(buffer[0], filelen[0], 1, fileptr[0]); // Read in the entire file
  fread(buffer[1], filelen[1], 1, fileptr[1]); // Read in the entire file
  fclose(fileptr[0]); // Close the file
  fclose(fileptr[1]); // Close the file

  std::vector<int> evals(c.num_wires);
  //Read input the right way!
  for (int i = 0; i < c.num_inp_wires; ++i) {
    evals[i] = GetBitReversed(i, buffer[0]);
  }

  //Evaluate the AND gates
  for (int i = 0; i < c.num_gates; ++i) {
    Gate g = c.gates[i];
    if (g.type == NOT) {
      evals[g.out_wire] = !evals[g.left_wire];
    } else if (g.type == XOR) {
      evals[g.out_wire] = evals[g.left_wire] ^ evals[g.right_wire];
    } else if (g.type == AND) {
      evals[g.out_wire] = evals[g.left_wire] & evals[g.right_wire];
    }
  }

  uint8_t* res = new uint8_t[16]();
  for (int i = 0; i < c.num_out_wires; ++i) {
    if (evals[c.num_wires - c.num_out_wires + i]) {
      SetBitReversed(i, 1, res);
    } else {
      SetBitReversed(i, 0, res);
    }
  }

  ASSERT_TRUE(std::equal(res, res + BITS_TO_BYTES(c.num_out_wires), buffer[1]));
}