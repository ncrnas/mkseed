project (test_findseed)

cmake_minimum_required (VERSION 2.8.2)

include_directories( ${gtest_SOURCE_DIR}/include
        ${gtest_SOURCE_DIR})

add_executable(findseed_test basic_check.cpp)

target_link_libraries(findseed_test gtest gtest_main)
