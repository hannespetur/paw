set(stations_test_files
  test_all_of.cpp
  test_any_of.cpp
  test_count_if.cpp
  test_count.cpp
  test_fill.cpp
  test_for_each.cpp
  test_internal.cpp
  test_none_of.cpp
  test_partition_iterator.cpp
  test_sort.cpp
  test_split.cpp
)

add_executable(test_pawstation ${stations_test_files})

# Link to threading library (pthreads on linux)
target_link_libraries (test_pawstation shared catch)
target_include_directories(test_pawparser PRIVATE ../include)

add_test(NAME test_pawstation COMMAND test_pawstation)
