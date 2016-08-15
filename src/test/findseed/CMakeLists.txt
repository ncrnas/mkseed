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
        ${gtest_SOURCE_DIR} ${MKS_HEADERS})

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})

# Add executable applications
add_executable (findseed_test basic_check.cpp)

# Add dependencies found by find_package (SeqAn)
target_link_libraries (findseed_test gtest gtest_main)
target_link_libraries (findseed_test msd ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

