cmake_minimum_required (VERSION 3.2)

add_executable(ex1_basic_usage ex1_basic_usage.cpp)

add_custom_target(examples COMMAND ex1_basic_usage) # Use 'make check' to get verbose test results
