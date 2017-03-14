#ifndef TINY_GARBLING_HALFGATES_H_
#define TINY_GARBLING_HALFGATES_H_

#include "tiny-util/typedefs.h"

//Convenience class for storing pointers
class HalfGates {
public:
  //The garbled tables
  uint8_t* T_E;
  uint8_t* T_G;

  //The output key
  uint8_t* out_key;

  //The solderings
  uint8_t* S_L;
  uint8_t* S_R;
  uint8_t* S_O;
};

#endif /* TINY_GARBLING_HALFGATES_H_ */