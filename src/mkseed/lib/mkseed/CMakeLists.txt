cmake_minimum_required (VERSION 2.8.2)

set(MKS_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/mkseed_core.cpp
        PARENT_SCOPE)

set(MKS_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE)