#include "commit/commit-scheme-rec.h"

CommitReceiver::CommitReceiver(Params& params, uint8_t seeds[], uint8_t choices[]) : CommitScheme(params), seeds(seeds), choices(choices) {

  //Initialize the matrix variables used for the matrix consisting of postulated values in BatchDecommit
  row_dim_values = AES_BITS;
  row_dim_values_bytes = AES_BYTES;
  transpose_matrix_values_size = row_dim_values_bytes * col_dim;
}

bool CommitReceiver::Commit() {
  ExpandAndTranspose();
  CheckbitCorrection();

  return ConsistencyCheck();
}

void CommitReceiver::ExpandAndTranspose() {

  //We first initialize num_blocks + 1 blocks to hold the share values. The +1 is due to our matrix transposition not being in-place, so we need a "scratch-pad" block. It is initially the 0 block that is used for this.
  for (int j = 0; j < num_blocks + 1; ++j) {
    matrices.emplace_back(std::make_unique<uint8_t[]>(transpose_matrix_size));
  }

  //Increment counter appropriately based on exec_id and expand seeds using PRNG
  PRNG rnd;
  int increase_counter_pr_block = CEIL_DIVIDE(col_dim_bytes, AES_BYTES * PIPELINES);
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    //First we increment the seeds to the current block.
    uint128_t seeds0 = uint8_tTOuint128_t(seeds + i * CSEC_BYTES);
    seeds0 += increase_counter_pr_block * num_blocks * params.exec_id;

    // Set seed
    rnd.SetSeed((uint8_t*) &seeds0);

    //Fill up the i'th row of all blocks. Notice the 0'th block is not filled up.
    for (int j = 0; j < num_blocks; ++j) {
      rnd.GenRnd(matrices[j + 1].get() + i * col_dim_bytes, col_dim_bytes);
    }
  }

  // Transpose the j+1 block and store in the j'th block. Results in block 0,...,num_blocks-1 contain the transposed matrices. Then point expanded data into the share arrays for easy access
  for (int j = 0; j < num_blocks; ++j) {
    transpose_320_128(matrices[j + 1].get(), matrices[j].get(), col_blocks);
    for (int i = 0; i < col_dim; ++i) {
      if (j * col_dim + i < num_commits_produced) { //last block might not be filled up
        commit_shares.emplace_back(matrices[j].get() + i * row_dim_bytes);
      }
    }
  }
}

void CommitReceiver::CheckbitCorrection() {
  //Receive correction values
  std::unique_ptr<uint8_t[]> checkbit_corrections_buf(std::make_unique<uint8_t[]>(num_commits_produced * BCH_BYTES));
  params.chan.ReceiveBlocking(checkbit_corrections_buf.get(), num_commits_produced * BCH_BYTES);

  //Run over all commitments and apply the correction to the checkbits. We do this using XOR and AND to do this efficiency. The correction is only applied if we hold the 1-share, so we use byte-wise ANDing to "select" the correction bits in each byte.
  for (int j = 0; j < num_commits_produced; ++j) {
    for (int p = 0; p < BCH_BYTES; ++p) {
      commit_shares[j][CSEC_BYTES + p] ^= (checkbit_corrections_buf[j * BCH_BYTES + p] & REVERSE_BYTE_ORDER[choices[CSEC_BYTES + p]]);
    }
  }
}

bool CommitReceiver::ConsistencyCheck() {
  //Sample and send the consistency check challenge element alpha.
  uint8_t alpha_seed[CSEC_BYTES];
  params.rnd.GenRnd(alpha_seed, CSEC_BYTES);
  params.chan.Send(alpha_seed, CSEC_BYTES);


  //Below is the naive approach for computing the linear combinations. Can be used to compare to the gf method
  // // ///////////////////////TEST//////////////////////////////////
  // auto t = GET_TIME();
  // uint8_t* V = new uint8_t[BITS_TO_BYTES(params.num_commits * CSEC)];
  // params.rnd.GenRnd(V, BITS_TO_BYTES(params.num_commits * CSEC));
  // uint8_t* linear_combination = new uint8_t[CODEWORD_BYTES*CSEC];
  // for (int l = 0; l < params.num_commits; l++) {
  //   for (int i = 0; i < CSEC; i++) {
  //     //Lookup in V to determine which indices go in each linear combination
  //     if (GetBit(l * CSEC + i, V)) {
  //       XOR_CodeWords(linear_combination +i*CODEWORD_BYTES, commit_shares[l]);//Builds them into the last num_cons_checks indices as these are only supposed to be used for this check anyway. After this they are "sacrificed" and therefore not reusable.
  //     }
  //   }
  // }
  // delete[] V;
  // delete[] linear_combination;
  // auto t_end = GET_TIME();
  // // ///////////////////////TEST//////////////////////////////////
  // auto t_new = GET_TIME();

  //Setup all registers for calculation the linear combinations. Will end up with 2*SSEC linear combinations.
  uint8_t final_result[CODEWORD_BYTES * 2 * SSEC];

  //res_tmp is twice as large as we do not do degree reduction until the very end, so we need to accumulate a larger intermediate value.
  __m128i res_tmp[2][CODEWORD_BITS];
  __m128i res_total[CODEWORD_BITS];
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    res_tmp[0][i] = _mm_setzero_si128();
    res_tmp[1][i] = _mm_setzero_si128();
    res_total[i] = _mm_setzero_si128();
  }

  __m128i val;
  __m128i val_result[2];

  //Need two temporary matrices for transposing each block of commitments which are added the the temporary results res_tmp.
  std::unique_ptr<uint8_t[]> matrices_tmp(std::make_unique<uint8_t[]>(2 * transpose_matrix_size));

  //Load the initial challenge element
  __m128i alpha = _mm_lddqu_si128((__m128i *) alpha_seed);
  for (int j = 0; j < num_blocks; ++j) {

    //Copy j'th block into temporary matrix and transpose. Result is in matrices_tmp.get().
    std::copy(matrices[j].get(), matrices[j].get() + transpose_matrix_size, matrices_tmp.get() + transpose_matrix_size);
    transpose_128_320(matrices_tmp.get() + transpose_matrix_size, matrices_tmp.get(), col_blocks);
    for (int l = 0; l < col_blocks; ++l) {
      for (int i = 0; i < CODEWORD_BITS; ++i) {
        //Pads transposed matrix with 0s if we are in the last block and it is not filled up. Needed to not destroy the linear combinations by reading garbage.
        if (j * col_dim + l * AES_BITS > num_commits_produced - AES_BITS) {
          int diff = j * col_dim + l * AES_BITS - (num_commits_produced  - AES_BITS);
          for (int p = 0; p < diff; ++p) {
            SetBitReversed(AES_BITS - diff + p, 0, matrices_tmp.get() + l * AES_BYTES + i * col_dim_bytes);
          }
        }
        //Load current row into val. If we are in one of the last AES_BITS commitments we directly add this to the final result as blinding. Else we multiply by alpha^(i+1) and store it in res_tmp.
        val = _mm_lddqu_si128((__m128i*) (matrices_tmp.get() + l * AES_BYTES + i * col_dim_bytes));

        if (j * col_dim + l * AES_BITS < num_commits_produced - AES_BITS) {
          //The actual commitments are multiplied with alpha
          mul128_karatsuba(val, alpha, &val_result[0], &val_result[1]);

          //Accumulate the val_result into res_tmp
          res_tmp[0][i] = _mm_xor_si128(res_tmp[0][i], val_result[0]);
          res_tmp[1][i] = _mm_xor_si128(res_tmp[1][i], val_result[1]);
        } else {
          //The AES_BITS blinding one-time commitments are added directly to res_total
          res_total[i] = val;
        }
      }
      //When done with one col_block we square the challenge element alpha. There are 8 col_blocks within each block
      gfmul128_no_refl(alpha, alpha, &alpha);
    }
  }
  //For comparing the time using the gf approach and above commented out naive approach.
  // auto t_new_end = GET_TIME();
  // PRINT_TIME(t_end, t, "naive");
  // PRINT_TIME(t_new_end, t_new, "gf");

  //mask is used to select the first 2*SSEC linear combinations from res_total and store in final_result. Needed as we actually produce AES_BITS linear combinations due to convenience. However we only send and verify 2*SSEC of these.
  uint8_t mask[CSEC_BYTES] = {0};
  std::fill(mask, mask + 2 * SSEC_BYTES, 0xFF);
  __m128i store_mask = _mm_lddqu_si128((__m128i*) mask);

  //Reduce and store the resulting linear combinations
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    gfred128_no_refl(res_tmp[0][i], res_tmp[1][i], &res_tmp[0][i]);
    res_total[i] = _mm_xor_si128(res_total[i], res_tmp[0][i]);

    //Finally move the 2*SSEC first linear combinations into final_result
    _mm_maskmoveu_si128(res_total[i], store_mask, (char*) (final_result + i * 2 * SSEC_BYTES));
  }

  //Receive the decommitments from CommitSnd that will be compared to the computed shares in final_result.
  uint8_t decommit_shares0[2 * CODEWORD_BITS * 2 * SSEC_BYTES];
  uint8_t* decommit_shares1 = decommit_shares0 + CODEWORD_BITS * 2 * SSEC_BYTES;
  params.chan.ReceiveBlocking(decommit_shares0, 2 * CODEWORD_BITS * 2 * SSEC_BYTES);

  if (!VerifyTransposedDecommits(decommit_shares0, decommit_shares1, final_result, 2 * SSEC)) {
    return false;
  }

  return true;

}

bool CommitReceiver::BatchDecommit(uint8_t computed_shares[], int num_values, uint8_t values[]) {

  //Sample and send the consistency check challenge element alpha.
  uint8_t alpha_seed[CSEC_BYTES];
  params.rnd.GenRnd(alpha_seed, CSEC_BYTES);
  params.chan.Send(alpha_seed, CSEC_BYTES);

  //Setup all registers for calculation the linear combinations. Will end up with SSEC linear combinations. The final_values_result will hold the linear combinations of the postulated values.
  uint8_t final_result[CODEWORD_BYTES * SSEC];
  uint8_t final_values_result[CSEC_BYTES * SSEC];

  //res_tmp is twice as large as we do not do degree reduction until the very end, so we need to accumulate a larger intermediate value.
  __m128i res_tmp[2][CODEWORD_BITS];
  __m128i res_values_tmp[2][AES_BITS];
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    res_tmp[0][i] = _mm_setzero_si128();
    res_tmp[1][i] = _mm_setzero_si128();
    if (i < AES_BITS) {
      res_values_tmp[0][i] = _mm_setzero_si128();
      res_values_tmp[1][i] = _mm_setzero_si128();
    }
  }

  __m128i val;
  __m128i val_result[2];

  //Used for transposing the commitment shares
  std::unique_ptr<uint8_t[]> matrices_tmp(std::make_unique<uint8_t[]>(2 * transpose_matrix_size));

  //Used for transposing the values
  CBitVector matrix_buffer_values;
  matrix_buffer_values.AttachBuf((uint8_t*) malloc(transpose_matrix_values_size), transpose_matrix_values_size);

  //Load the initial challenge element
  __m128i alpha = _mm_lddqu_si128((__m128i *) alpha_seed);

  //Compute number of check_blocks needed in total for num_values
  int num_check_blocks = CEIL_DIVIDE(num_values, col_dim);

  //For each check_block we load the shares in column-major order and then transpose to get to row-major order so we can address AES_BITS values entry-wise at a time.

  for (int j = 0; j < num_check_blocks; ++j) {
    //Load block
    for (int i = 0; i < col_dim; ++i) {

      //First we check if we are done, ie we are in the last block and it is not entirely filled up.
      int num_check_index = j * col_dim + i;
      if (num_check_index < num_values) {
        //We copy the correct share into our current block matrix
        std::copy(computed_shares + num_check_index * CODEWORD_BYTES, computed_shares + num_check_index * CODEWORD_BYTES + CODEWORD_BYTES, matrices_tmp.get() + transpose_matrix_size + i * row_dim_bytes);

        //We copy the current value into our current block matrix
        std::copy(values + num_check_index * CSEC_BYTES, values + num_check_index * CSEC_BYTES + CSEC_BYTES, matrix_buffer_values.GetArr() + i * row_dim_values_bytes);

      } else {
        //This pads the last block with 0 rows. Can be optimized by only doing this once. However need to calculate the remaining number of "rows" to zero out then.
        std::fill(matrices_tmp.get() + transpose_matrix_size + i * row_dim_bytes, matrices_tmp.get() + transpose_matrix_size + i * row_dim_bytes + CODEWORD_BYTES, 0);

        std::fill(matrix_buffer_values.GetArr() + i * row_dim_values_bytes, matrix_buffer_values.GetArr() + i * row_dim_values_bytes + CSEC_BYTES, 0);
      }
    }

    //Transpose block
    transpose_128_320(matrices_tmp.get() + transpose_matrix_size, matrices_tmp.get(), col_blocks);

    //The values matrix can be transposed using 8*128x128 matrix transpose.
    matrix_buffer_values.EklundhBitTranspose(col_dim, row_dim_values);

    //Compute on block. Processes the block matrices in the same way as for the consistency check with the modification that there is no blinding values
    for (int l = 0; l < col_blocks; ++l) {
      for (int i = 0; i < CODEWORD_BITS; ++i) {
        val = _mm_lddqu_si128((__m128i*) (matrices_tmp.get() + l * AES_BYTES + i * col_dim_bytes));
        mul128_karatsuba(val, alpha, &val_result[0], &val_result[1]);

        res_tmp[0][i] = _mm_xor_si128(res_tmp[0][i], val_result[0]);
        res_tmp[1][i] = _mm_xor_si128(res_tmp[1][i], val_result[1]);

        //Also, only process the values for the first AES_BITS bits as these are only this long
        if (i < AES_BITS) {
          val = _mm_lddqu_si128((__m128i*) (matrix_buffer_values.GetArr() + l * AES_BYTES + i * col_dim_bytes));
          mul128_karatsuba(val, alpha, &val_result[0], &val_result[1]);
          res_values_tmp[0][i] = _mm_xor_si128(res_values_tmp[0][i], val_result[0]);
          res_values_tmp[1][i] = _mm_xor_si128(res_values_tmp[1][i], val_result[1]);
        }
      }
      //When done with one col_block we square the challenge element alpha. There are 8 col_blocks within each check_block
      gfmul128_no_refl(alpha, alpha, &alpha);
    }
  }

  //Need to manually delete this matrix
  matrix_buffer_values.delCBitVector();

  //mask is used to select the first SSEC linear combinations from res_total and store in final_result. Needed as we actually produce AES_BITS linear combinations due to convenience. However we only send and verify 2*SSEC of these.
  uint8_t mask[CSEC_BYTES] = {0};
  std::fill(mask, mask + SSEC_BYTES, 0xFF);
  __m128i store_mask = _mm_lddqu_si128((__m128i*) mask);

  //Reduce and store the resulting linear combinations
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    gfred128_no_refl(res_tmp[0][i], res_tmp[1][i], &res_tmp[0][i]);

    //Finally move the SSEC first linear combinations of the shares into final_result
    _mm_maskmoveu_si128(res_tmp[0][i], store_mask, (char*) (final_result + i * SSEC_BYTES));

    if (i < AES_BITS) {
      gfred128_no_refl(res_values_tmp[0][i], res_values_tmp[1][i], &res_values_tmp[0][i]);

      //Finally move the SSEC first linear combinations of the values into final_values_result
      _mm_maskmoveu_si128(res_values_tmp[0][i], store_mask, (char*) (final_values_result + i * SSEC_BYTES));
    }
  }

  //Receive the decommitments from CommitSnd that will be compared to the computed shares in final_result.
  uint8_t decommit_shares0[2 * CODEWORD_BYTES * SSEC];
  uint8_t* decommit_shares1 = decommit_shares0 + CODEWORD_BYTES * SSEC;
  params.chan.ReceiveBlocking(decommit_shares0, 2 * CODEWORD_BYTES * SSEC);

  //Verify the received shares are valid
  if (!VerifyTransposedDecommits(decommit_shares0, decommit_shares1, final_result, SSEC)) {
    return false;
  }

  //If the decommits check out, processed and verify that these were actually linear combinations of the postulated values as well. If everything is ok we are sure that the num_values postulated values in values[] are indeed the committed values and we return true.
  for (int i = 0; i < SSEC; ++i) {
    XOR_128(decommit_shares0 + i * AES_BYTES, decommit_shares1 + i * AES_BYTES);

    if (!std::equal(final_values_result + i * AES_BYTES, final_values_result + i * AES_BYTES + AES_BYTES, decommit_shares0 + i * AES_BYTES)) {
      return false;
    }
  }

  return true;
}

//This only implements a single call to chosen commit and it is very coupled to the TinyLEGO application. If needed it should be pretty easy to add support for multiple invocations and generalize it.
void CommitReceiver::ChosenCommit(int num_values) {
  num_chosen_commits = num_values;

  chosen_commit_values = std::move(std::unique_ptr<uint8_t[]>(std::make_unique<uint8_t[]>(num_values * CSEC_BYTES)));

  params.chan.ReceiveBlocking(chosen_commit_values.get(), num_values * CSEC_BYTES);
}

//Implements a single chosen decommit call. The positions of idxs are into the chosen_commit_values array and decides which of these are used for de-masking. Should not be hard to generalize to multiple calls.
bool CommitReceiver::ChosenDecommit(uint8_t computed_shares[], uint8_t decommit_values_res[], std::vector<uint64_t> idxs, int num_values) {

  std::unique_ptr<uint8_t[]> decommit_shares0(std::make_unique<uint8_t[]>(num_values * (CODEWORD_BYTES + CSEC_BYTES)));
  uint8_t* decommit_shares1 = decommit_shares0.get() + num_values * CODEWORD_BYTES;

  params.chan.ReceiveBlocking(decommit_shares0.get(), num_values * (CODEWORD_BYTES + CSEC_BYTES));

  if (!VerifyDecommits(decommit_shares0.get(), decommit_shares1, computed_shares, decommit_values_res, choices, code.get(), num_values)) {
    return false;
  }

  //The decommitted values are stored in decommit_values_res and we now unmask them using the previously sent correction values in chosen_commit_values
  for (int i = 0; i < num_values; ++i) {
    XOR_128(decommit_values_res + i * CSEC_BYTES, chosen_commit_values.get() + idxs[i] * CSEC_BYTES);
  }

  return true;
}

bool CommitReceiver::VerifyTransposedDecommits(uint8_t decommit_shares0[], uint8_t decommit_shares1[], uint8_t computed_shares[], int num_values) {

  if (num_values > AES_BITS) {
    throw std::runtime_error("Unsupported number of values");
  }

  std::unique_ptr<uint8_t[]> matrix0(std::make_unique<uint8_t[]>(2 * transpose_matrix_size));
  uint8_t* matrix1 = matrix0.get() + transpose_matrix_size;

  //Used to select only the first num_values_bytes of each share. Gives us a way of only reading the first SSEC or 2*SSEC linear combinations of an otherwise 128-bit register and zero'ing out the remaining unused bits.
  int num_values_bytes = BITS_TO_BYTES(num_values);
  uint8_t mask[CSEC_BYTES] = {0};
  std::fill(mask, mask + num_values_bytes, 0xFF);
  __m128i store_mask = _mm_lddqu_si128((__m128i*) mask);

  //Read all shares in row-major order and compare choice bits row-wise. We check all num_values in parallel this way. Notice we select only the first num_values bits of each share using store_mask and bit-wise AND.
  __m128i share, share0, share1, tmp;
  for (int i = 0; i < CODEWORD_BITS; ++i) {
    share = _mm_lddqu_si128((__m128i*) (computed_shares + i * num_values_bytes));
    share = _mm_and_si128(share, store_mask);

    share0 = _mm_lddqu_si128((__m128i*) (decommit_shares0 + i * num_values_bytes));
    share0 = _mm_and_si128(share0, store_mask);

    share1 = _mm_lddqu_si128((__m128i*) (decommit_shares1 + i * num_values_bytes));
    share1 = _mm_and_si128(share1, store_mask);

    tmp = _mm_xor_si128(share0, share1);
    _mm_storeu_si128((__m128i*) (matrix0.get() + i * col_dim_single_bytes), tmp);

    if (GetBit(i, choices)) {
      if (!compare128(share, share1)) {
        return false;

      }
    } else {
      if (!compare128(share, share0)) {
        return false;
      }
    }
  }

  //At this point we know that the decommitments match the computed_shares in each position. Next we need to verify that the linear combinations are themselves codewords. Only if both requirements hold we are sure the decommitment is valid.

  //Transpose the decommitted values so we can access them column-wise and thus compute the check-bits.
  transpose_320_128(matrix0.get(), matrix1);

  uint8_t tmp_checkbits[BCH_BYTES];
  for (int i = 0; i < num_values; ++i) {

    //Reset the tmp_checkbits in every check
    std::fill(tmp_checkbits, tmp_checkbits + BCH_BYTES, 0);

    //Compute the actual codeword of the first CSEC bits and store in tmp_checkbits.
    code->Encode(matrix1 + i * row_dim_bytes, tmp_checkbits);

    //Ensure that tmp_checkbits is actually equal to the checkbits of the column. This ensures that entire column is a codeword.
    if (!std::equal(tmp_checkbits, tmp_checkbits + BCH_BYTES, matrix1 + CSEC_BYTES + i * row_dim_bytes)) {
      std::cout << "Abort! Linear combination " << i << " is not a codeword" << std::endl;
      return false; //Not a codeword!
    }
  }

  return true; //All checks passed!
}