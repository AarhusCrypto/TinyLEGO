#ifndef TINY_COMMIT_ECC_H_
#define TINY_COMMIT_ECC_H_

#include "tiny-util/util.h"
extern "C" {
  #include "commit/bch.h"
}

//This class acts as a wrapper for invoking the BCH code
class ECC {
public:
  ECC();
  void Encode(uint8_t data[], uint8_t checkbits[]);
  
  std::unique_ptr<struct bch_control> bch_control;
};

#endif /* TINY_COMMIT_ECC_H_ */