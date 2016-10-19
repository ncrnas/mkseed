# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /test/findseed/CMakeLists.txt
#
# CMakeLists.txt file for find_seed
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (test_findseed)
message (STATUS "Configuring test/find_seed")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search dependencies
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../mkseed/lib/mkseed msd)

# Add include directories
include_directories (${SEQAN_INCLUDE_DIRS} ${gtest_SOURCE_DIR}/include
        ${gtest_SOURCE_DIR} ${MKS_HEADERS}
        ${CMAKE_CURRENT_LIST_DIR}/../include)

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})
if (MSVC)
    add_definitions (/D TEST_DATA_PATH=${CMAKE_CURRENT_LIST_DIR}/../data/)
else()
    add_definitions (-DTEST_DATA_PATH=${CMAKE_CURRENT_LIST_DIR}/../data/)
endif()

# Add executable applications
add_executable (findseed_test basic_check.cpp read_fasta.cpp
        u3_001.cpp u3_002.cpp u3_003.cpp u3_004.cpp u3_005.cpp u3_006.cpp u3_007.cpp u3_008.cpp
        u3_009.cpp u3_010.cpp u3_011.cpp
        mi_001.cpp)

# Add dependencies found by find_package (SeqAn)
target_link_libraries (findseed_test gtest gtest_main)
target_link_libraries (findseed_test msd ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Define helper functions and macros used by Google Test.
include(${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest/cmake/internal_utils.cmake)
config_compiler_and_linker()  # Defined in internal_utils.cmake.
