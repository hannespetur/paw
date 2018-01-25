cmake_minimum_required (VERSION 2.8.11)

add_executable(ex1 ex1_basic_usage.cpp)
add_executable(ex2 ex2_simple.cpp)
add_executable(ex3 ex3_subcommand.cpp)

if (Boost_FILESYSTEM_FOUND AND Boost_SYSTEM_FOUND)
  add_executable(ex4 ex4_boost_filesystem.cpp)
  target_link_libraries(ex4 ${Boost_LIBRARIES})
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(ex4 shared)
endif()


foreach(index RANGE 1 3)
  target_link_libraries("ex${index}" shared)
endforeach()
