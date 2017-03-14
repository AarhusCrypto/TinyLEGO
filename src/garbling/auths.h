#ifndef TINY_GARBLING_AUTHS_H_
#define TINY_GARBLING_AUTHS_H_

#include "tiny-util/typedefs.h"

//Convenience class for storing pointers
class Auths {
public:
  //The hashed values
  uint8_t* H_0;
  uint8_t* H_1;

  //Soldering
  uint8_t* S_A;
};

#endif /* TINY_GARBLING_AUTHS_H_ */