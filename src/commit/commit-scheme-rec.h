#ifndef TINY_COMMIT_COMMITSCHEME_REC_H_
#define TINY_COMMIT_COMMITSCHEME_REC_H_

#include "commit/commit-scheme.h"

class CommitReceiver : public CommitScheme {
public:
  CommitReceiver(Params& params, uint8_t seeds[], uint8_t choices[]);

  //Committing
  bool Commit();
  void ExpandAndTranspose();
  void CheckbitCorrection();
  bool ConsistencyCheck();

  //Chosen Commit/Decommit.
  void ChosenCommit(int num_values);
  bool ChosenDecommit(uint8_t computed_shares[], uint8_t decommit_values_res[], std::vector<uint64_t> idxs, int num_values);

  //BatchDecommit
  bool BatchDecommit(uint8_t computed_shares[], int num_values, uint8_t values[]);

  //Verify
  bool VerifyTransposedDecommits(uint8_t decommit_shares0[], uint8_t decommit_shares1[], uint8_t computed_shares[], int num_values);

  //Size of the matrices used for transposing the postulated opening values
  int row_dim_values;
  int row_dim_values_bytes;
  int transpose_matrix_values_size;

  //Holds the actual data
  std::vector<std::unique_ptr<uint8_t[]>> matrices;

  //Convenience pointers
  uint8_t* seeds;
  uint8_t* choices;
  std::vector<uint8_t*> commit_shares;

  //Chosen commits data
  int num_chosen_commits;
  std::unique_ptr<uint8_t[]> chosen_commit_values;

};

//This function is static inline for efficiency reasons as it's used in the online phase of TinyLEGO (and ChosenDecommit, but the reason it's static is due to the online phase part).
static inline bool VerifyDecommits(uint8_t decommit_share0[], uint8_t decommit_share1[], uint8_t computed_shares[], uint8_t res_values[], uint8_t choices[], ECC* code, int num_values) {

  uint8_t c[BCH_BYTES] = {0};
  uint8_t c1[BCH_BYTES];

  for (int j = 0; j < num_values; ++j) {
    //Check value shares
    for (int i = 0; i < CSEC_BYTES; ++i) {
      if (((computed_shares + j * CODEWORD_BYTES)[i] ^ ((decommit_share1 + j * CSEC_BYTES)[i] & REVERSE_BYTE_ORDER[choices[i]]) ^ ((decommit_share0 + j * CODEWORD_BYTES)[i] & ~REVERSE_BYTE_ORDER[choices[i]])) != 0) {
        return false;
      }
    }

    //Construct c1 checkbit shares
    XOR_128(res_values + j * CSEC_BYTES, decommit_share0 + j * CODEWORD_BYTES, decommit_share1 + j * CSEC_BYTES);
    code->Encode(res_values + j * CSEC_BYTES, c);
    XOR_CheckBits(c1, c, (decommit_share0 + j * CODEWORD_BYTES) + CSEC_BYTES);
    std::fill(c, c + BCH_BYTES, 0);

    //Check checkbit shares
    for (int i = 0; i < BCH_BYTES; ++i) {
      if (((computed_shares + j * CODEWORD_BYTES)[CSEC_BYTES + i] ^ (c1[i] & REVERSE_BYTE_ORDER[choices[CSEC_BYTES + i]]) ^ ((decommit_share0 + j * CODEWORD_BYTES)[CSEC_BYTES + i] & ~REVERSE_BYTE_ORDER[choices[CSEC_BYTES + i]])) != 0) {
        return false;
      }
    }
  }

  return true; //All checks passed!
};

#endif /* TINY_COMMIT_COMMITSCHEME_REC_H_ */