#ifndef TINY_DOT_EXT_H_
#define TINY_DOT_EXT_H_


#include "OTExtension/ot/alsz-ot-ext-rec.h"
#include "OTExtension/ot/alsz-ot-ext-snd.h"
#include "OTExtension/ot/ot-ext-snd.h"
#include "OTExtension/ot/ot-ext-rec.h"
#include "OTExtension/ot/xormasking.h"

#include "tiny/params.h"

class ALSZDOTExt {
public:
  ALSZDOTExt(Params& params);

  Params& params;
  Network net;

  int bit_length_outer;
  uint32_t num_seed_OT;
  int bit_length_inner;
  uint32_t num_check_OT;

  int num_snd_vals;
  field_type m_eFType;

  uint32_t num_OT_threads;
  snd_ot_flavor s_type;
  rec_ot_flavor r_type;

  bool m_bUseMinEntCorAssumption;
  std::unique_ptr<MaskingFunction> m_fMaskFct;

protected:
  void PrivacyAmplification(uint8_t priv_amp_matrix[], int rows_bytes, int columns_bits, int num_vecs, uint8_t base_inner[], uint8_t base_outer[]);
  void GeneratePrivAmpMatrix(uint8_t priv_amp_seed[], uint8_t priv_amp_matrix[], int size);
};

#endif /* TINY_DOT_EXT_H_ */