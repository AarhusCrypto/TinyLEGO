#include "commit/commit-scheme.h"

CommitScheme::CommitScheme(Params& params) : params(params), code(std::make_unique<ECC>()) {
  
  row_dim = PAD_TO_MULTIPLE(CODEWORD_BITS, COMMIT_HALF_BLOCK_SIZE); //Makes row_dim 320 which is needed for our matrix transposition code.
  row_dim_bytes = BITS_TO_BYTES(row_dim);
  col_blocks = COMMIT_NUM_TRANSPOSE_BLOCKS; //Chosen for efficiency. This gave the best result.
  col_dim_single = AES_BITS; //We process AES_BITS columns at a time with our PCLMULQDQ code.
  col_dim_single_bytes = BITS_TO_BYTES(col_dim_single);  
  col_dim = col_blocks * col_dim_single;
  col_dim_bytes = BITS_TO_BYTES(col_dim);

  num_commits_produced = params.num_commits + AES_BITS; //We produce AES_BITS extra commitments that we use for blinding. However only the first 2*SSEC (Consistency Check) and SSEC (BatchDecommit) are actually sent over the network and checked.

  num_blocks = CEIL_DIVIDE(num_commits_produced, col_dim);
  transpose_matrix_size = BITS_TO_BYTES(row_dim * col_dim);
}