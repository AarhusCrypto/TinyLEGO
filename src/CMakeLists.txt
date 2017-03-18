set(PARAMS_SRCS tiny/params.cpp)
add_library(PARAMS ${PARAMS_SRCS})
target_link_libraries(PARAMS SPLIT_COMMIT)

set(CIRCUIT_SRCS circuit/circuit-parser.cpp)
add_library(CIRCUIT ${CIRCUIT_SRCS})

set(GARBLING_SRCS garbling/garbling-handler.cpp)
add_library(GARBLING ${GARBLING_SRCS})

set(TINY_SRCS tiny/tiny-evaluator.cpp tiny/tiny-constructor.cpp tiny/tiny.cpp)
add_library(TINY ${TINY_SRCS})
target_link_libraries(TINY PARAMS GARBLING CIRCUIT)

add_executable(Tinyconst mains/tiny-const-main.cpp)
target_link_libraries(Tinyconst TINY)

add_executable(Tinyeval mains/tiny-eval-main.cpp)
target_link_libraries(Tinyeval TINY)