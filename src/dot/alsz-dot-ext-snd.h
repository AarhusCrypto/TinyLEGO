#ifndef TINY_DOT_EXT_SND_H_
#define TINY_DOT_EXT_SND_H_

#include "dot/alsz-dot-ext.h"

class ALSZDOTExtSnd : public ALSZDOTExt {
public:
  ALSZDOTExtSnd(Params& params, bool set_lsb_delta);

  void InitOTSender();

  void Send();
  void PrivacyAmplification(uint8_t base_inner[], uint8_t delta_inner[]);

  ALSZOTExtSnd sender;

  std::unique_ptr<uint8_t[]> base_outer;
  std::unique_ptr<uint8_t[]> delta_outer;
  bool set_lsb_delta;
};

#endif /* TINY_DOT_EXT_SND_H_ */