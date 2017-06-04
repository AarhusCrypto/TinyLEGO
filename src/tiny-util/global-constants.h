#ifndef TINY_UTIL_GLOBAL_CONSTANTS_H_
#define TINY_UTIL_GLOBAL_CONSTANTS_H_

#include "SplitCommit/src/util/global-constants.h"

//Behavior
// #define DEBUG_SOLDERINGS_INP_BUCKETS

static uint8_t tiny_constant_seeds[2][16] = {
  {0x43, 0x73, 0x98, 0x42, 0x70, 0x13, 0x38, 0x78, 0xAB, 0x45, 0x78, 0xFF, 0xEA, 0xD3, 0xFF, 0x03},
  {0x43, 0x73, 0x98, 0x41, 0x70, 0x12, 0x38, 0x78, 0x43, 0x73, 0x98, 0x41, 0x66, 0x19, 0xAA, 0xFE}};

#define CSEC 128
#define CSEC_BYTES 16
#define SSEC 40
#define SSEC_BYTES 5

static uint8_t global_aes_key[CSEC_BYTES] = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};

//ThreadPool
#define NUM_IO_THREADS 1
#define TP_MUL_FACTOR 2

//q, alpha, beta, p_a, p_g
#define bucket_param_table_size 45
static int bucket_param_table[bucket_param_table_size][5] = {
  {96900000, 3, 3, 3, 3},
  {78150000, 3, 3, 2, 3},
  {64150000, 3, 4, 7, 8},
  {44450000, 3, 4, 8, 7},
  {31700000, 3, 4, 7, 7},
  {25300000, 3, 4, 6, 7},
  {22100000, 3, 4, 5, 7},
  {17500000, 3, 4, 7, 6},
  {12450000, 3, 4, 6, 6},
  {9950000, 3, 4, 5, 6},
  {8650000, 3, 4, 4, 6},
  {6850000, 3, 4, 6, 5},
  {4850000, 3, 4, 5, 5},
  {3850000, 3, 4, 4, 5},
  {2650000, 3, 4, 5, 4},
  {1850000, 3, 4, 4, 4},
  {1450000, 3, 4, 3, 4},
  {1000000, 3, 4, 4, 3},
  {700000, 3, 4, 3, 3},
  {550000, 3, 4, 2, 3,},
  {350000, 3, 4, 3, 2},
  {221696, 3, 4, 2, 2},
  {200000, 4, 5, 5, 4},
  {150000, 4, 5, 4, 4},
  {100000, 4, 5, 3, 4},
  {91081, 4, 5, 4, 3},
  {50000, 4, 5, 3, 3},
  {37460, 4, 5, 2, 3},
  {34000, 5, 6, 5, 4},
  {25000, 4, 5, 3, 2},
  {19500, 5, 6, 3, 4},
  {16500, 4, 5, 2, 2},
  {10000, 5, 6, 3, 3},
  {8000, 5, 6, 2, 3},
  {6928, 6, 7, 3, 4},
  {5500, 5, 6, 3, 2},
  {5000, 6, 7, 4, 3},
  {3500, 5, 6, 2, 2},
  {3000, 6, 7, 2, 3},
  {2500, 7, 8, 4, 3},
  {2000, 6, 7, 3, 2},
  {1500, 6, 7, 2, 2},
  {1000, 7, 8, 3, 2},
  {500, 8, 9, 3, 2}
  {7, 40, 41, 1, 1}
};

#define bucket_param_table_online_size 12
static int bucket_param_table_online[bucket_param_table_online_size][5] = {
  {96900000, 3, 3, 3, 3},
  {78150000, 3, 3, 2, 3},
  {73400000, 4, 3, 8, 3},
  {46350000, 3, 3, 3, 2},
  {27200000, 3, 3, 2, 2},
  {21350000, 4, 3, 8, 2},
  {20650000, 4, 3, 7, 2},
  {11050000, 3, 3, 2, 1},
  {6950000, 4, 3, 8, 1},
  {4650000, 4, 3, 7, 1},
  {4100000, 4, 3, 6, 1},
  {4000000, 4, 3, 5, 1}
};
#endif /* TINY_UTIL_GLOBAL_CONSTANTS_H_ */
