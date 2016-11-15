#ifndef TINY_CIRCUIT_CIRCUIT_H_
#define TINY_CIRCUIT_CIRCUIT_H_

#include "util/typedefs.h"

enum GATE {
  AND = 0,
  OR = 1,
  XOR = 2,
  NAND = 3,
  NXOR = 4,
  NOR = 5,
  RXOR = 6,
  LXOR = 7,
  NOT = 8
};

class Gate {
public:
  uint32_t left_wire;
  uint32_t right_wire;
  uint32_t out_wire;
  enum GATE type;
};

class Circuit {
public:
  std::vector<Gate> gates;
  uint32_t num_wires;
  uint32_t num_const_inp_wires;
  uint32_t num_eval_inp_wires;
  uint32_t num_inp_wires;
  uint32_t num_out_wires;
  uint32_t num_and_gates;
  uint32_t num_gates;
};

Circuit read_text_circuit(const char* circuit_file);
Circuit ParseCircuit(char* data);

#endif /* TINY_CIRCUIT_CIRCUIT_H_ */