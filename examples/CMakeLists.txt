cmake_minimum_required (VERSION 2.8)

add_library(impl OBJECT implementions.cpp)

add_executable(ex1 $<TARGET_OBJECTS:impl> ex1_basic_usage.cpp)
add_executable(ex2 $<TARGET_OBJECTS:impl> ex2_simple.cpp)
add_executable(ex3 $<TARGET_OBJECTS:impl> ex3_subcommand.cpp)
