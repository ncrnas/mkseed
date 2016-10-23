cmake_minimum_required (VERSION 2.8.2)
project (mkseed)

# Add include paths
set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/seqan/extras/include
                                              ${CMAKE_SOURCE_DIR}/seqan/core/include
                                              ${CMAKE_SOURCE_DIR})

# Add apps subdirectory
add_subdirectory (apps)

