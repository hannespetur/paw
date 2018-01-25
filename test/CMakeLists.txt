cmake_minimum_required (VERSION 2.8.11)

add_library(catch OBJECT main.cpp)

add_executable(test_pawparser $<TARGET_OBJECTS:catch> test_pawparser.cpp)
target_link_libraries(test_pawparser shared)

add_test(NAME pawparser_catch_tests COMMAND $<TARGET_FILE:test_pawparser>)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose) # Use 'make check' to get verbose test results
