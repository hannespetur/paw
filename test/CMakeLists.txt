cmake_minimum_required (VERSION 3.2)

add_executable(test_pawparser test_pawparser.cpp)

enable_testing()
add_test(NAME pawparser_catch_tests COMMAND test_pawparser)
