#ifndef TINY_DOT_EXT_REC_H_
#define TINY_DOT_EXT_REC_H_

#include "dot/alsz-dot-ext.h"

class ALSZDOTExtRec : public ALSZDOTExt {
public:
  ALSZDOTExtRec(Params& params);
  void InitOTReceiver();

  void Receive();
  void PrivacyAmplification(uint8_t response_inner[]);

  ALSZOTExtRec receiver;
  std::unique_ptr<uint8_t[]> choices_outer;
  std::unique_ptr<uint8_t[]> response_outer;
};

#endif /* TINY_DOT_EXT_REC_H_ */