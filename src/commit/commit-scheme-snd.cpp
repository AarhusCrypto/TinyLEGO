#include "commit/commit-scheme-snd.h"

CommitSender::CommitSender(Params& params, uint8_t seeds0[], uint8_t seeds1[]) : CommitScheme(params), seeds0(seeds0), seeds1(seeds1) {
}

void CommitSender::Commit() {
  ExpandAndTranspose();
  CheckbitCorrection();
  ConsistencyCheck();
}

void CommitSender::ExpandAndTranspose() {

  //We first initialize num_blocks + 1 blocks to hold the share values. The +1 is due to our matrix transposition not being in-place, so we need a "scratch-pad" block. It is initially the 0 block that is used for this.
  for (int j = 0; j < num_blocks + 1; ++j) {
    matrices0.emplace_back(std::make_unique<uint8_t[]>(transpose_matrix_size));
    matrices1.emplace_back(std::make_unique<uint8_t[]>(transpose_matrix_size));
  }

  //Increment counter appropriately based on exec_id and expand seeds using PRNGs
  PRNG rnd0;
  PRNG rnd1;
  int increase_counter_pr_block = CEIL_DIVIDE(col_dim_bytes, AES_BYTES * PIPELINES);
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    //First we increment the seeds to the current block.
    uint128_t seed0 = uint8_tTOuint128_t(seeds0 + i * CSEC_BYTES);
    uint128_t seed1 = uint8_tTOuint128_t(seeds1 + i * CSEC_BYTES);
    seed0 += increase_counter_pr_block * num_blocks * params.exec_id;
    seed1 += increase_counter_pr_block * num_blocks * params.exec_id;

    // Set seed
    rnd0.SetSeed((uint8_t*)&seed0);
    rnd1.SetSeed((uint8_t*)&seed1);

    //Fill up the i'th row of all blocks. Notice the 0'th block is not filled up.
    for (int j = 0; j < num_blocks; ++j) {
      rnd0.GenRnd(matrices0[j + 1].get() + i * col_dim_bytes, col_dim_bytes);
      rnd1.GenRnd(matrices1[j + 1].get() + i * col_dim_bytes, col_dim_bytes);
    }
  }

  // Transpose the j+1 block and store in the j'th block. Results in block 0,...,num_blocks-1 contain the transposed matrices. Then point expanded data into the shares arrays for easy access
  for (int j = 0; j < num_blocks; ++j) {
    transpose_320_128(matrices0[j + 1].get(), matrices0[j].get(), col_blocks);
    transpose_320_128(matrices1[j + 1].get(), matrices1[j].get(), col_blocks);
    for (int i = 0; i < col_dim; ++i) {
      if (j * col_dim + i < num_commits_produced) { //last block might not be filled up
        commit_shares0.emplace_back(matrices0[j].get() + i * row_dim_bytes);
        commit_shares1.emplace_back(matrices1[j].get() + i * row_dim_bytes);
      }
    }
  }
}

void CommitSender::CheckbitCorrection() {
  //Buffers
  std::unique_ptr<uint8_t[]> checkbit_corrections_buf(std::make_unique<uint8_t[]>(num_commits_produced * BCH_BYTES));
  uint8_t values_buffer[CODEWORD_BYTES];

  //Compute all the corrections and send to receiver
  for (int j = 0; j < num_commits_produced; ++j) {
    XOR_CodeWords(values_buffer, commit_shares0[j], commit_shares1[j]);
    code->Encode(values_buffer, checkbit_corrections_buf.get() + j * BCH_BYTES);

    XOR_CheckBits(commit_shares1[j] + CSEC_BYTES, commit_shares0[j] + CSEC_BYTES, checkbit_corrections_buf.get() + j * BCH_BYTES);
    XOR_CheckBits(checkbit_corrections_buf.get() + j * BCH_BYTES, values_buffer + CSEC_BYTES);
  }
  //Needs to be SendBlocking, else code hangs. Cannot explain why as everywhere else Send works fine.
  params.chan.SendBlocking(checkbit_corrections_buf.get(), num_commits_produced * BCH_BYTES);
}

void CommitSender::ConsistencyCheck() {

  //Setup all registers for calculation the linear combinations. Will end up with 2*SSEC_BITS linear combinations. Everything here is 2x larger than in commit-scheme-rec.cpp since we need to compute linear combinations of both shares.
  uint8_t final_result0[2 * CODEWORD_BYTES * 2 * SSEC];
  uint8_t* final_result1 = final_result0 + CODEWORD_BYTES * 2 * SSEC;

  //res_tmps is twice as large as we do not do degree reduction until the very end, so we need to accumulate a larger intermediate value.
  __m128i res_tmps[4][CODEWORD_BITS];
  __m128i res_totals[2][CODEWORD_BITS];
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    res_tmps[0][i] = _mm_setzero_si128();
    res_tmps[1][i] = _mm_setzero_si128();
    res_tmps[2][i] = _mm_setzero_si128();
    res_tmps[3][i] = _mm_setzero_si128();
    res_totals[0][i] = _mm_setzero_si128();
    res_totals[1][i] = _mm_setzero_si128();
  }

  __m128i vals[2];
  __m128i vals_result[4];

  //Need four temporary matrices for transposing each block of commitments which are added the the temporary results res_tmp. Each share needs two matrices.
  std::unique_ptr<uint8_t[]> matrices_tmp0(std::make_unique<uint8_t[]>(4 * transpose_matrix_size));
  uint8_t* matrices_tmp1 = matrices_tmp0.get() + 2 * transpose_matrix_size;

  //Receive challenge seed from receiver and load initial challenge alpha
  uint8_t alpha_seed[CSEC_BYTES];
  params.chan.ReceiveBlocking(alpha_seed, CSEC_BYTES);
  __m128i alpha = _mm_lddqu_si128((__m128i *) alpha_seed);

  for (int j = 0; j < num_blocks; ++j) {
    //Copy j'th block into temporary matrix and transpose. Result is in matrices_tmp0.get() and matrices_tmp1.
    std::copy(matrices0[j].get(), matrices0[j].get() + transpose_matrix_size, matrices_tmp0.get() + transpose_matrix_size);
    std::copy(matrices1[j].get(), matrices1[j].get() + transpose_matrix_size, matrices_tmp1 + transpose_matrix_size);
    transpose_128_320(matrices_tmp0.get() + transpose_matrix_size, matrices_tmp0.get(), col_blocks);
    transpose_128_320(matrices_tmp1 + transpose_matrix_size, matrices_tmp1, col_blocks);

    for (int l = 0; l < col_blocks; ++l) {
      for (int i = 0; i < CODEWORD_BITS; ++i) {
        //Pads matrix with 0s if we are in the last block and it is not filled up. Needed to not destroy the linear combinations by reading garbage.
        if (j * col_dim + l * AES_BITS > num_commits_produced - AES_BITS) {
          int diff = j * col_dim + l * AES_BITS - (num_commits_produced  - AES_BITS);
          for (int p = 0; p < diff; ++p) {
            SetBitReversed(AES_BITS - diff + p, 0, matrices_tmp0.get() + l * AES_BYTES + i * col_dim_bytes);
            SetBitReversed(AES_BITS - diff + p, 0, matrices_tmp1 + l * AES_BYTES + i * col_dim_bytes);
          }
        }
        //Load current row into val. If we are in one of the last AES_BITS commitments we directly add this to the final result as blinding. Else we multiply by alpha^(i+1) and store it in res_tmp.
        vals[0] = _mm_lddqu_si128((__m128i*) (matrices_tmp0.get() + l * AES_BYTES + i * col_dim_bytes));
        vals[1] = _mm_lddqu_si128((__m128i*) (matrices_tmp1 + l * AES_BYTES + i * col_dim_bytes));

        if (j * col_dim + l * AES_BITS < num_commits_produced - AES_BITS) {
          //The actual commitments are multiplied with alpha
          mul128_karatsuba(vals[0], alpha, &vals_result[0], &vals_result[1]);
          mul128_karatsuba(vals[1], alpha, &vals_result[2], &vals_result[3]);

          //Accumulate the vals_result into res_tmps
          res_tmps[0][i] = _mm_xor_si128(res_tmps[0][i], vals_result[0]);
          res_tmps[1][i] = _mm_xor_si128(res_tmps[1][i], vals_result[1]);
          res_tmps[2][i] = _mm_xor_si128(res_tmps[2][i], vals_result[2]);
          res_tmps[3][i] = _mm_xor_si128(res_tmps[3][i], vals_result[3]);
        } else {
          //The AES_BITS blinding one-time commitments are added directly to res_totals
          res_totals[0][i] = vals[0];
          res_totals[1][i] = vals[1];
        }
      }
      //When done with one col_block we square the challenge element alpha. There are 8 col_blocks within each block
      gfmul128_no_refl(alpha, alpha, &alpha);
    }
  }

  //mask is used to select the first 2*SSEC linear combinations from res_totals and store in final_result0 and final_results1. Needed as we actually produce AES_BITS linear combinations due to convenience. However we only send and verify 2*SSEC of these.
  uint8_t mask[CSEC_BYTES] = {0};
  std::fill(mask, mask + 2 * SSEC_BYTES, 0xFF);
  __m128i store_mask = _mm_lddqu_si128((__m128i*) mask);

  //Reduce and store the resulting linear combinations
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    gfred128_no_refl(res_tmps[0][i], res_tmps[1][i], &res_tmps[0][i]);
    gfred128_no_refl(res_tmps[2][i], res_tmps[3][i], &res_tmps[2][i]);
    res_totals[0][i] = _mm_xor_si128(res_totals[0][i], res_tmps[0][i]);
    res_totals[1][i] = _mm_xor_si128(res_totals[1][i], res_tmps[2][i]);

    //Finally move the 2*SSEC first linear combinations into final_result0 and final_result1
    _mm_maskmoveu_si128(res_totals[0][i], store_mask, (char*) (final_result0 + i * 2 * SSEC_BYTES));
    _mm_maskmoveu_si128(res_totals[1][i], store_mask, (char*) (final_result1 + i * 2 * SSEC_BYTES));
  }
  //Send the resulting 2*SSEC decommitments
  params.chan.Send(final_result0, 2 * CODEWORD_BYTES * 2 * SSEC);
}

void CommitSender::BatchDecommit(uint8_t decommit_shares0[], uint8_t decommit_shares1[], int num_values) {

  //Setup all registers for calculation the linear combinations. Will end up with SSEC_BITS linear combinations. Everything here is 2x larger than in commit-scheme-rec.cpp since we need to compute linear combinations of both shares.
  uint8_t final_result0[2 * CODEWORD_BYTES * SSEC];
  uint8_t* final_result1 = final_result0 + CODEWORD_BYTES * SSEC;

  //res_tmps is twice as large as we do not do degree reduction until the very end, so we need to accumulate a larger intermediate value. As opposed to ConsistencyCheck, we do not need res_totals here as we do not use any blinding values in BatchDecommit.
  __m128i res_tmps[4][CODEWORD_BITS];
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    res_tmps[0][i] = _mm_setzero_si128();
    res_tmps[1][i] = _mm_setzero_si128();
    res_tmps[2][i] = _mm_setzero_si128();
    res_tmps[3][i] = _mm_setzero_si128();
  }

  __m128i vals[2];
  __m128i vals_result[4];

  //Need four temporary matrices for transposing each block of commitments which are added the the temporary results res_tmp. Each share needs two matrices.
  std::unique_ptr<uint8_t[]> matrices_tmp0(std::make_unique<uint8_t[]>(4 * transpose_matrix_size));
  uint8_t* matrices_tmp1 = matrices_tmp0.get() + 2 * transpose_matrix_size;

  //Receive challenge seed from receiver and load initial challenge alpha
  uint8_t alpha_seed[CSEC_BYTES];
  params.chan.ReceiveBlocking(alpha_seed, CSEC_BYTES);
  __m128i alpha = _mm_lddqu_si128((__m128i *) alpha_seed);

  //Compute number of check_blocks needed in total for num_values
  int num_check_blocks = CEIL_DIVIDE(num_values, col_dim);

  //For each check_block we load the decommitments in column-major order and then transpose to get to row-major order so we can address AES_BITS values entry-wise at a time.

  for (int j = 0; j < num_check_blocks; ++j) {
    //Load block
    for (int i = 0; i < col_dim; ++i) {
      int num_check_index = j * col_dim + i;
      if (num_check_index < num_values) {
        std::copy(decommit_shares0 + num_check_index * CODEWORD_BYTES, decommit_shares0 + num_check_index * CODEWORD_BYTES + CODEWORD_BYTES, matrices_tmp0.get() + transpose_matrix_size + i * row_dim_bytes);
        std::copy(decommit_shares1 + num_check_index * CODEWORD_BYTES, decommit_shares1 + num_check_index * CODEWORD_BYTES + CODEWORD_BYTES, matrices_tmp1 + transpose_matrix_size + i * row_dim_bytes);
      }
      else { //this pads the last block with 0 rows
        std::fill(matrices_tmp0.get() + transpose_matrix_size + i * row_dim_bytes, matrices_tmp0.get() + transpose_matrix_size + i * row_dim_bytes + CODEWORD_BYTES, 0);
        std::fill(matrices_tmp1 + transpose_matrix_size + i * row_dim_bytes, matrices_tmp1 + transpose_matrix_size + i * row_dim_bytes + CODEWORD_BYTES, 0);
      }
    }

    //Transpose block
    transpose_128_320(matrices_tmp0.get() + transpose_matrix_size, matrices_tmp0.get(), col_blocks);
    transpose_128_320(matrices_tmp1 + transpose_matrix_size, matrices_tmp1, col_blocks);

    //Compute on block. Processes the block matrices in the same way as for the consistency check with the modification that there is no blinding values
    for (int l = 0; l < col_blocks; ++l) {
      for (int i = 0; i < CODEWORD_BITS; ++i) {
        vals[0] = _mm_lddqu_si128((__m128i*) (matrices_tmp0.get() + l * AES_BYTES + i * col_dim_bytes));
        vals[1] = _mm_lddqu_si128((__m128i*) (matrices_tmp1 + l * AES_BYTES + i * col_dim_bytes));
        mul128_karatsuba(vals[0], alpha, &vals_result[0], &vals_result[1]);
        mul128_karatsuba(vals[1], alpha, &vals_result[2], &vals_result[3]);

        res_tmps[0][i] = _mm_xor_si128(res_tmps[0][i], vals_result[0]);
        res_tmps[1][i] = _mm_xor_si128(res_tmps[1][i], vals_result[1]);
        res_tmps[2][i] = _mm_xor_si128(res_tmps[2][i], vals_result[2]);
        res_tmps[3][i] = _mm_xor_si128(res_tmps[3][i], vals_result[3]);
      }
      //When done with one col_block we square the challenge element alpha. There are 8 col_blocks within each check_block
      gfmul128_no_refl(alpha, alpha, &alpha);
    }
  }

  //mask is used to select the first SSEC linear combinations from res_totals and store in final_result0 and final_results1. Needed as we actually produce AES_BITS linear combinations due to convenience. However we only send and verify SSEC of these.
  uint8_t mask[CSEC_BYTES] = {0};
  std::fill(mask, mask + SSEC_BYTES, 0xFF);
  __m128i store_mask = _mm_lddqu_si128((__m128i*) mask);

  //Reduce and store the resulting linear combinations
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    gfred128_no_refl(res_tmps[0][i], res_tmps[1][i], &res_tmps[0][i]);
    gfred128_no_refl(res_tmps[2][i], res_tmps[3][i], &res_tmps[2][i]);

     //Finally move the SSEC first linear combinations into final_result0 and final_result1
    _mm_maskmoveu_si128(res_tmps[0][i], store_mask, (char*) (final_result0 + i * SSEC_BYTES));
    _mm_maskmoveu_si128(res_tmps[2][i], store_mask, (char*) (final_result1 + i * SSEC_BYTES));
  }

  //Send the resulting SSEC decommitments
  params.chan.Send(final_result0, 2 * CODEWORD_BYTES * SSEC);
}

//Implements a single chosen decommit using random commitment positions in idxs. Should not be hard to generalize to multiple calls.
void CommitSender::ChosenCommit(uint8_t values[], std::vector<uint64_t> idxs, int num_values) {

  std::unique_ptr<uint8_t[]> tmp_values(std::make_unique<uint8_t[]>(num_values * CSEC_BYTES));

  //Mask the values using the predescribed commitments with index indxs[i].
  for (int i = 0; i < num_values; ++i) {
    XOR_128(tmp_values.get() + i * CSEC_BYTES, values + i * CSEC_BYTES, commit_shares0[idxs[i]]);
    XOR_128(tmp_values.get() + i * CSEC_BYTES, commit_shares1[idxs[i]]);
  }

  params.chan.Send(tmp_values.get(), num_values * CSEC_BYTES);
}