#ifndef TINY_COMMIT_COMMITSCHEME_SND_H_
#define TINY_COMMIT_COMMITSCHEME_SND_H_

#include "commit/commit-scheme.h"

class CommitSender : public CommitScheme {
public:
  CommitSender(Params& params, uint8_t seeds0[], uint8_t seeds1[]);

  void Commit();
  void ChosenCommit(uint8_t values[], std::vector<uint64_t> idxs, int num_values);
  void BatchDecommit(uint8_t decommit_shares0[], uint8_t decommit_shares1[], int num_values);
    
  //Holds the actual data
  std::vector<std::unique_ptr<uint8_t[]>> matrices0;
  std::vector<std::unique_ptr<uint8_t[]>> matrices1;

  //Convenience pointers
  uint8_t* seeds0;
  uint8_t* seeds1;
  std::vector<uint8_t*> commit_shares0;
  std::vector<uint8_t*> commit_shares1;

private:
  void ExpandAndTranspose();
  void CheckbitCorrection();
  void ConsistencyCheck();
};
#endif /* TINY_COMMIT_COMMITSCHEME_SND_H_ */