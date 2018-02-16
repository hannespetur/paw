cmake_minimum_required (VERSION 2.8.11)

add_executable(test_align test_align.cpp)
target_link_libraries(test_align shared ${PAW_LINK_LIBRARIES})

add_test(NAME test_pawalign COMMAND $<TARGET_FILE:test_align>)
