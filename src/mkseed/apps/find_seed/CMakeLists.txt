cmake_minimum_required (VERSION 2.8.2)
project (mkseed_apps_find_seed)

message (STATUS "Configuring mikan/apps/find_seed")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search dependencies
set (SEQAN_FIND_DEPENDENCIES "NONE")
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------
option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." OFF)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/mkseed msd)

# Add include directories
include_directories (${SEQAN_INCLUDE_DIRS}
        ${MKS_HEADERS})

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})

# Add executable applications
add_executable (find_seed find_seed.cpp)

# Add dependencies found by find_package (SeqAn)
add_library(msd ${MKS_SOURCE})

target_link_libraries (find_seed msd ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Define helper functions and macros used by Google Test.
include(${CMAKE_CURRENT_SOURCE_DIR}/../../../test/lib/googletest/cmake/internal_utils.cmake)
config_compiler_and_linker()  # Defined in internal_utils.cmake.

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Install find_seed in ${PREFIX}/bin directory
install (TARGETS find_seed
         DESTINATION bin)

install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/../../../../LICENSE
         ${CMAKE_CURRENT_SOURCE_DIR}/../../../../README.md
         DESTINATION share/doc/mkseed)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

set (CPACK_PACKAGE_NAME "find_seed")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "find_seed: identify 6mer miRNA target sites on given seqeuences")
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Takaya Saito <takaya.saito@outlook.com>")
set (CPACK_PACKAGE_VENDOR "Takaya Saito <takaya.saito@outlook.com>")
set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/../../../../LICENSE)
set (CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/../../../../README.md)
set (CPACK_PACKAGE_VERSION "${MKSEED_VERSION}")
set (CPACK_PACKAGE_VERSION_MAJOR "${MKSEED_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${MKSEED_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${MKSEED_VERSION_PATCH}")

include (CPack)
