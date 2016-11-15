#include "circuit/circuit.h"

//Parse the gate description given a char array of the description file.
Circuit ParseCircuit(char raw_circuit[]) {
  Circuit circuit;
  raw_circuit = strchr(raw_circuit, ' ') + 1; //we dont need num_gates
  circuit.num_wires = (uint32_t) atoi(raw_circuit);
  raw_circuit = strchr(raw_circuit,  '\n') + 1; //Skip to next line

  circuit.num_const_inp_wires = (uint32_t) atoi(raw_circuit);
  raw_circuit = strchr(raw_circuit, ' ') + 1;
  circuit.num_eval_inp_wires = (uint32_t) atoi(raw_circuit);
  raw_circuit = strchr(raw_circuit, ' ') + 1;
  circuit.num_out_wires = (uint32_t) atoi(raw_circuit);
  circuit.num_inp_wires = circuit.num_const_inp_wires + circuit.num_eval_inp_wires;

  raw_circuit = strchr(raw_circuit,  '\n') + 1; //Skip to next line
  raw_circuit = strchr(raw_circuit,  '\n') + 1; //Skip to next line
  int curr_gate_num = 0;
  circuit.num_and_gates = 0;
  uint32_t num_inputs, left_wire_idx, right_wire_idx, out_wire_idx;
  char type[3];

  while (*raw_circuit != EOF) {
    if (*raw_circuit == '\n') {
      raw_circuit = strchr(raw_circuit,  '\n') + 1;
      continue;
    }
    num_inputs = (uint32_t) atoi(raw_circuit);
    raw_circuit = strchr(raw_circuit,  ' ') + 1;
    raw_circuit = strchr(raw_circuit,  ' ') + 1; //We skip num_output wires as they all have 1.

    if (num_inputs == 1) {
      left_wire_idx = (uint32_t) atoi(raw_circuit);
      raw_circuit = strchr(raw_circuit,  ' ') + 1;
      out_wire_idx = (uint32_t) atoi(raw_circuit);
      raw_circuit = strchr(raw_circuit,  ' ') + 1;
      raw_circuit = strchr(raw_circuit,  '\n') + 1;
      circuit.gates.emplace_back(Gate());
      circuit.gates[curr_gate_num].type = NOT;
      circuit.gates[curr_gate_num].left_wire = left_wire_idx;
      circuit.gates[curr_gate_num].out_wire = out_wire_idx;
      ++curr_gate_num;
    } else {
      left_wire_idx = (uint32_t) atoi(raw_circuit);
      raw_circuit = strchr(raw_circuit,  ' ') + 1;
      right_wire_idx = (uint32_t) atoi(raw_circuit);
      raw_circuit = strchr(raw_circuit,  ' ') + 1;
      out_wire_idx = (uint32_t) atoi(raw_circuit);

      raw_circuit = strchr(raw_circuit,  ' ') + 1;

      memcpy(type, raw_circuit, 3 * sizeof(char));
      std::string type_string(type);
      raw_circuit = strchr(raw_circuit,  '\n') + 1;
      if (type_string.find("XOR") != std::string::npos) {
        circuit.gates.emplace_back(Gate());
        circuit.gates[curr_gate_num].type = XOR;
        circuit.gates[curr_gate_num].left_wire = left_wire_idx;
        circuit.gates[curr_gate_num].right_wire = right_wire_idx;
        circuit.gates[curr_gate_num].out_wire = out_wire_idx;
        ++curr_gate_num;
      } else if (type_string.find("AND") != std::string::npos) {
        circuit.gates.emplace_back(Gate());
        circuit.gates[curr_gate_num].type = AND;
        circuit.gates[curr_gate_num].left_wire = left_wire_idx;
        circuit.gates[curr_gate_num].right_wire = right_wire_idx;
        circuit.gates[curr_gate_num].out_wire = out_wire_idx;

        ++curr_gate_num;
        ++circuit.num_and_gates;
      }
    }
  }

  //Add identity AND-gates to all output wires. This is to simplify TinyLEGO evaluation as now it is easy to identify which output wire needs to be leaked for the evalutator to decode the output.
  int out_start = circuit.num_wires - circuit.num_out_wires;
  int org_num_wires = circuit.num_wires;
  for (int i = 0; i < circuit.num_out_wires; ++i) {
    circuit.gates.emplace_back(Gate());
    circuit.gates[curr_gate_num].type = AND;
    circuit.gates[curr_gate_num].left_wire = out_start + i;
    circuit.gates[curr_gate_num].right_wire = out_start + i;
    circuit.gates[curr_gate_num].out_wire = org_num_wires + i;
    
    ++curr_gate_num;
    ++circuit.num_and_gates;
    ++circuit.num_wires;
  }

  circuit.num_gates = curr_gate_num;

  return circuit;
}

//Reads circuit in textual format. Writes byte length of text file to file_size.
Circuit read_text_circuit(const char* circuit_file) {
  FILE* file;
  size_t file_size;
  file = fopen(circuit_file, "r");
  if (file == NULL) {
    printf("ERROR: Could not open text circuit: %s\n", circuit_file);
    exit(EXIT_FAILURE);
  }
  fseek(file, 0, SEEK_END);
  file_size = ftell(file);
  rewind(file);

  std::unique_ptr<char[]> data(new char[file_size + 1]);
  size_t size = fread(data.get(), 1, file_size, file);
  if (size != file_size) {
    printf("ERROR while loading file: %s\n", circuit_file);
    exit(EXIT_FAILURE);
  }
  data[file_size] = EOF;
  if (ferror(file)) {
    printf("ERROR: fread() error\n");
    exit(EXIT_FAILURE);
  }
  fclose(file);
  Circuit circuit = ParseCircuit(data.get());

  return circuit;

}