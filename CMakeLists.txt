cmake_minimum_required (VERSION 2.8.3)
project (TINY C CXX)

#Use -DCMAKE_BUILD_TYPE=RELEASE to invoke
IF(${CMAKE_BUILD_TYPE} MATCHES "RELEASE")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build/release/")
set(CMAKE_VERBOSE_MAKEFILE FALSE)

#Defaults to debug mode
ELSE(${CMAKE_BUILD_TYPE} MATCHES "RELEASE")
set(${CMAKE_BUILD_TYPE} "DEBUG")
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build/debug/")
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "RELEASE")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build/obj/")

#Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -g3 -O0 -Wall -pg")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG -g3 -O0 -Wall -pg")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DRELEASE")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -DRELEASE")

#Find libraries

set(PTHREAD_COMP_FLAG -pthread)
FIND_PACKAGE(OpenSSL REQUIRED)
FIND_PACKAGE(FindGMP HINTS cmake/)
FIND_PACKAGE(FindZMQ HINTS cmake/)

#Required for large mallocs on some platforms
set(ENV{MALLOC_PERMIT_INSANE_REQUESTS} "SOME VALUE")

#set constant compiler flags
set(CMAKE_C_FLAGS "-march=native ${PTHREAD_COMP_FLAG}")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++14")

#Include headers
include_directories(${GMP_INCLUDE_DIR})

#Build sublibraries
add_subdirectory(libs/SplitCommit)

#Include headers
include_directories(${CMAKE_CURRENT_LIST_DIR}/libs/SplitCommit/libs/libOTe)
include_directories(${CMAKE_CURRENT_LIST_DIR}/libs/SplitCommit/src/)
include_directories(${CMAKE_CURRENT_LIST_DIR}/libs/SplitCommit/libs/)

include_directories(${OPENSSL_INCLUDE_DIR})
include_directories(${ZMQ_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test)



#Build sublibraries
##################################MIRACL#######################################
set(MIRACL_OTX_SRCS
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrcore.c
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrarth0.c
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrarth1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrarth2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mralloc.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrsmall.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrio1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrio2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrgcd.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrjack.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrxgcd.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrarth3.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrbits.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrrand.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrprime.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrcrt.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrscrt.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrmonty.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrpower.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrsroot.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrcurve.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrfast.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrshs.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrshs256.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrshs512.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mraes.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrgcm.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrlucas.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrzzn2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrzzn2b.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrzzn3.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrzzn4.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrecn2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrstrong.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrbrick.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrebrick.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrec2m.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrgf2m.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrflash.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrfrnd.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrdouble.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrround.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrbuild.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrflsh1.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrpi.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrflsh2.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrflsh3.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrflsh4.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/mrmuldv.c 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/big.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/zzn.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/ecn.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/ec2.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/crt.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/flash.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/polymod.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/libs/Miracl/poly.cpp)
add_library(MIRACL_OTX ${MIRACL_OTX_SRCS})
##################################MIRACL_OTX#######################################

#################################ALSZ OTX#####################################
set(OTX_UTIL_SRCS libs/OTExtension/util/cbitvector.cpp libs/OTExtension/util/parse_options.cpp)
add_library(OTX_UTIL ${OTX_UTIL_SRCS})
set(CRYPTO_SRCS libs/OTExtension/util/crypto/crypto.cpp libs/OTExtension/util/crypto/ecc-pk-crypto.cpp libs/OTExtension/util/crypto/gmp-pk-crypto.cpp libs/OTExtension/util/crypto/intrin_sequential_enc8.cpp libs/OTExtension/util/crypto/TedKrovetzAesNiWrapperC.cpp)
add_library(OTX_CRYPTO ${CRYPTO_SRCS})
target_link_libraries(OTX_CRYPTO MIRACL_OTX ${OPENSSL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

set(OTX_SRCS libs/OTExtension/ot/alsz-ot-ext-rec.cpp libs/OTExtension/ot/alsz-ot-ext-snd.cpp libs/OTExtension/ot/ot-ext-rec.cpp libs/OTExtension/ot/iknp-ot-ext-rec.cpp libs/OTExtension/ot/iknp-ot-ext-snd.cpp libs/OTExtension/ot/kk-ot-ext-rec.cpp libs/OTExtension/ot/kk-ot-ext-snd.cpp libs/OTExtension/ot/nnob-ot-ext-rec.cpp libs/OTExtension/ot/nnob-ot-ext-snd.cpp libs/OTExtension/ot/ot-ext-snd.cpp libs/OTExtension/ot/ot-ext.cpp libs/OTExtension/ot/naor-pinkas.cpp libs/OTExtension/ot/pvwddh.cpp libs/OTExtension/ot/simpleot.cpp)
add_library(OTX ${OTX_SRCS})
target_link_libraries(OTX OTX_UTIL OTX_CRYPTO)
#################################ALSZ OTX#####################################

#add sources
add_subdirectory(src)

#add tests
add_subdirectory(test)