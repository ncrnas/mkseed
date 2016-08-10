# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /mikan/CMakeLists.txt
#
# CMakeLists.txt file for mikan
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (mikan)

# Add include paths
set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../seqan/extras/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/../seqan/core/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Setup the library modules from core and extras
seqan_setup_library (core)

# Add apps subdirectory
add_subdirectory (apps)

