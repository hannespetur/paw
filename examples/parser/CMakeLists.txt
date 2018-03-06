add_executable(parser_ex1 ex1_basic_usage.cpp)
add_executable(parser_ex2 ex2_simple.cpp)
add_executable(parser_ex3 ex3_subcommand.cpp)

if (Boost_FILESYSTEM_FOUND AND Boost_SYSTEM_FOUND)
  add_executable(parser_ex4 ex4_boost_filesystem.cpp)
  target_link_libraries(parser_ex4 shared)
endif()

foreach(index RANGE 1 3)
  target_link_libraries("parser_ex${index}" shared)
endforeach()
