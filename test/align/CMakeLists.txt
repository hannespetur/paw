add_executable(test_align test_align.cpp)

target_link_libraries(test_align shared catch)
target_include_directories(test_align PRIVATE ../include)

add_test(NAME test_pawalign COMMAND $<TARGET_FILE:test_align>)
