#add_subdirectory(algorithm)
#add_subdirectory(helper)
#add_subdirectory(meta)
#add_subdirectory(memory)
#add_subdirectory(pack)
#add_subdirectory(range)

add_executable(test_boost_simd_api
    algorithm/test_all_of.cpp
    algorithm/test_any_of.cpp
    algorithm/test_copy.cpp
    algorithm/test_copy_n.cpp
    algorithm/test_count.cpp
    algorithm/test_count_if.cpp
    algorithm/test_equal.cpp
    algorithm/test_fill.cpp
    algorithm/test_find.cpp
    algorithm/test_find_if.cpp
    algorithm/test_find_if_not.cpp
    algorithm/test_generate.cpp
    algorithm/test_generate_n.cpp
    algorithm/test_iota.cpp
    algorithm/test_lexicographical_compare.cpp
    algorithm/test_max_element.cpp
    algorithm/test_max_val.cpp
    algorithm/test_min_element.cpp
    algorithm/test_min_val.cpp
    algorithm/test_none_of.cpp
    algorithm/test_not_all_of.cpp
    algorithm/test_reduce.cpp
    algorithm/test_replace.cpp
    algorithm/test_replace_if.cpp
    algorithm/test_transform.cpp
    algorithm/test_transform_reduce.cpp
    helper/test_logical.cpp
    helper/test_mask.cpp
    )

target_link_libraries(test_boost_simd_api shared catch)
target_include_directories(test_boost_simd_api PRIVATE ${PAW_INCLUDE_DIRS})
add_test(NAME test_boost_simd_api COMMAND $<TARGET_FILE:test_boost_simd_api>)
