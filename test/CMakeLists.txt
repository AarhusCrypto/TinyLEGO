# add_executable(TestCommitment test-commitment.cpp)
# target_link_libraries(TestCommitment COMMIT PARAMS DOT gtest_main gtest)

# add_executable(TestDOT test-dot.cpp)
# target_link_libraries(TestDOT PARAMS DOT gtest_main gtest)

# add_executable(TestDOTAndCommit test-dot-and-commit.cpp)
# target_link_libraries(TestDOTAndCommit COMMIT PARAMS DOT gtest_main gtest)

add_executable(TestParser test-circuit-parser.cpp)
target_link_libraries(TestParser CIRCUIT gtest_main gtest)

add_executable(TestTiny test-tiny.cpp)
target_link_libraries(TestTiny TINY gtest_main gtest)