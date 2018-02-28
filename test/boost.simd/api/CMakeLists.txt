#add_subdirectory(algorithm)
#add_subdirectory(helper)
#add_subdirectory(meta)
#add_subdirectory(memory)
#add_subdirectory(pack)
#add_subdirectory(range)

add_executable(test_boost_simd_algorithm
    algorithm.cpp
    algorithm/test_all_of.cpp
    function/test_abs.cpp
    )

target_link_libraries(test_boost_simd_algorithm shared ${PAW_LINK_LIBRARIES})
add_test(NAME test_boost_simd_algorithm COMMAND $<TARGET_FILE:test_boost_simd_algorithm>)
