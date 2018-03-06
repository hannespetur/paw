add_executable(test_pawparser test_pawparser.cpp)
target_link_libraries(test_pawparser shared catch)
target_include_directories(test_pawparser PRIVATE ../include)

add_test(NAME test_pawparser COMMAND $<TARGET_FILE:test_pawparser>)
