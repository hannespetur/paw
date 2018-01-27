cmake_minimum_required (VERSION 2.8.11)

add_executable(test_pawparser test_pawparser.cpp)
target_link_libraries(test_pawparser shared ${PAW_LINK_LIBRARIES})

add_test(NAME test_pawparser COMMAND $<TARGET_FILE:test_pawparser>)
