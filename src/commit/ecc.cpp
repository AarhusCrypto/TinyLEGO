#include "commit/ecc.h"

ECC::ECC() : bch_control((init_bch(CONFIG_BCH_CONST_M, CONFIG_BCH_CONST_T, 0))) {

}

//checkbits should be BCH_BYTES long and initialized to 0!
void ECC::Encode(uint8_t data[], uint8_t checkbits[]) {
  encode_bch(bch_control.get(), data, CSEC_BYTES, checkbits);
}