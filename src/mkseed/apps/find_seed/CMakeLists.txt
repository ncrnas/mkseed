# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /mikan/apps/find_seed/CMakeLists.txt
#
# CMakeLists.txt file for find_seed
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (mikan_apps_find_seed)
message (STATUS "Configuring mikan/apps/find_seed")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search dependencies
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})

# Add executable applications
add_executable (find_seed find_seed.cpp)

# Add dependencies found by find_package (SeqAn)
target_link_libraries (find_seed ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (find_seed)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install find_seed in ${PREFIX}/bin directory
install (TARGETS find_seed
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/find_seed for SeqAn release builds
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (find_seed)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:find_seed")
  set (CPACK_PACKAGE_NAME "find_seed")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "find_seed")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (find_seed "find_seed")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:find_seed")
