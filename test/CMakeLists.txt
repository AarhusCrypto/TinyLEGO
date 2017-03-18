add_executable(TestParser test-circuit-parser.cpp)
target_link_libraries(TestParser CIRCUIT gtest_main gtest)

add_executable(TestTiny test-tiny.cpp)
target_link_libraries(TestTiny TINY gtest_main gtest)