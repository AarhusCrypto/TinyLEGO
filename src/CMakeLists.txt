set(PARAMS_SRCS tiny/params.cpp)
add_library(PARAMS ${PARAMS_SRCS})
target_link_libraries(PARAMS SPLIT_COMMIT)

set(CIRCUIT_SRCS circuit/circuit-parser.cpp)
add_library(CIRCUIT ${CIRCUIT_SRCS})
target_link_libraries(CIRCUIT SPLIT_COMMIT)

set(GARBLING_SRCS garbling/garbling-handler.cpp)
add_library(GARBLING ${GARBLING_SRCS})
target_link_libraries(GARBLING SPLIT_COMMIT)

set(TINY_SRCS tiny/tiny-evaluator.cpp tiny/tiny-constructor.cpp tiny/tiny.cpp)
add_library(TINY ${TINY_SRCS})
target_link_libraries(TINY PARAMS GARBLING CIRCUIT)

add_executable(Tinyconst tiny-mains/tiny-const-main.cpp)
target_link_libraries(Tinyconst TINY)

add_executable(Tinyeval tiny-mains/tiny-eval-main.cpp)
target_link_libraries(Tinyeval TINY)