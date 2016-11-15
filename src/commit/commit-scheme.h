#ifndef TINY_COMMIT_COMMITSCHEME_H_
#define TINY_COMMIT_COMMITSCHEME_H_

#include "tiny/params.h"
#include "commit/ecc.h"

class CommitScheme {
public:
  CommitScheme(Params& params);

  Params& params;
  
  std::unique_ptr<ECC> code;

  uint64_t num_commits_produced;

  int row_dim;
  int row_dim_bytes;
  int col_dim_single;
  int col_dim_single_bytes;
  int col_blocks;
  int col_dim;
  int col_dim_bytes;
  int transpose_matrix_size;
  uint64_t num_blocks;
};

#endif /* TINY_COMMIT_COMMITSCHEME_H_ */