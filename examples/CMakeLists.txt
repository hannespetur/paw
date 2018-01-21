cmake_minimum_required (VERSION 2.8)

add_library(impl OBJECT implementions.cpp)

add_executable(ex1 $<TARGET_OBJECTS:impl> ex1_basic_usage.cpp)
add_executable(ex2 $<TARGET_OBJECTS:impl> ex2_simple.cpp)
add_executable(ex3 $<TARGET_OBJECTS:impl> ex3_subcommand.cpp)

if (Boost_FILESYSTEM_FOUND AND Boost_SYSTEM_FOUND)
  message (STATUS "Boost filesystem was found!")

  add_executable(ex4 $<TARGET_OBJECTS:impl> ex4_boost_filesystem.cpp)
  target_link_libraries(ex4 ${Boost_LIBRARIES})
  include_directories(${Boost_INCLUDE_DIRS})
endif()
