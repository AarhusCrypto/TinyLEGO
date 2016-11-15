set(NETWORK_SRCS util/network.cpp)
add_library(NETWORK ${NETWORK_SRCS})

set(CHANNEL_SRCS util/channel.cpp)
add_library(CHANNEL ${CHANNEL_SRCS})
target_link_libraries(CHANNEL ${ZMQ_LIBRARIES})

set(BCH_SRCS commit/bch.c)
add_library(BCH ${BCH_SRCS})

set(PRG_SRCS prg/random.cpp prg/aes-ni.cpp)
add_library(PRG ${PRG_SRCS})

set(PARAMS_SRCS tiny/params.cpp)
add_library(PARAMS ${PARAMS_SRCS})
target_link_libraries(PARAMS NETWORK OTX_CRYPTO PRG CHANNEL)

set(CIRCUIT_SRCS circuit/circuit-parser.cpp)
add_library(CIRCUIT ${CIRCUIT_SRCS})

set(DOT_SRCS dot/alsz-dot-ext-rec.cpp dot/alsz-dot-ext-snd.cpp dot/alsz-dot-ext.cpp)
add_library(DOT ${DOT_SRCS})
target_link_libraries(DOT NETWORK OTX)

set(COMMIT_SRCS commit/commit-scheme-rec.cpp commit/commit-scheme-snd.cpp commit/commit-scheme.cpp commit/ecc.cpp)
add_library(COMMIT ${COMMIT_SRCS})
target_link_libraries(COMMIT BCH PRG)

set(GARBLING_SRCS garbling/garbling-handler.cpp)
add_library(GARBLING ${GARBLING_SRCS})

set(TINY_SRCS tiny/tiny-evaluator.cpp tiny/tiny-constructor.cpp tiny/tiny.cpp)
add_library(TINY ${TINY_SRCS})
target_link_libraries(TINY COMMIT PARAMS DOT GARBLING CIRCUIT)

add_executable(Tinyconst mains/tiny-const-main.cpp)
target_link_libraries(Tinyconst TINY)

add_executable(Tinyeval mains/tiny-eval-main.cpp)
target_link_libraries(Tinyeval TINY)

add_executable(Commitsnd mains/commit-snd-main.cpp)
target_link_libraries(Commitsnd DOT COMMIT PARAMS)

add_executable(Commitrec mains/commit-rec-main.cpp)
target_link_libraries(Commitrec DOT COMMIT PARAMS)