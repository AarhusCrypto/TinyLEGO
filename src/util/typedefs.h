#ifndef TINY_UTIL_TYPEDEFS_H_
#define TINY_UTIL_TYPEDEFS_H_


// #include <mmintrin.h>  //MMX
// #include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
// #include <pmmintrin.h> //SSE3
// #include <tmmintrin.h> //SSSE3
// #include <smmintrin.h> //SSE4.1
// #include <nmmintrin.h> //SSE4.2
// #include <ammintrin.h> //SSE4A
#include <wmmintrin.h> //AES
#include <immintrin.h> //AVX
// #include <zmmintrin.h> //AVX512

#include <CTPL-master/ctpl_stl.h>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <climits>
#include <vector>
#include <thread>
#include <chrono>
#include <algorithm>
#include <numeric>

typedef unsigned __int128 uint128_t;
typedef uint8_t octet;
typedef unsigned int uint;

#endif /* TINY_UTIL_TYPEDEFS_H_ */