cmake_minimum_required(VERSION 3.20)
project(DMet)

set(CMAKE_CXX_STANDARD 14)

#make the dmet library
set(SOURCES
        src/DMet/DistribDistances.cpp
        src/DMet/PointDistances.cpp
        src/DMet/utils.cpp
        src/DMet/EqWidthBin.cpp
        )

include_directories(include/)
add_library(DMet SHARED ${SOURCES})


#include gmp headers and mpfr
include_directories("/usr/local/include/")
target_link_libraries(DMet PUBLIC "/usr/local/Cellar/gmp/6.2.1_1/lib/libgmp.a")
target_link_libraries(DMet PUBLIC "/usr/local/Cellar/mpfr/4.1.0/lib/libmpfr.a")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#Add test subdirectory
add_subdirectory(test)

option(BUILD-DOCS "BUILD DOCS" ON)
if(BUILD-DOCS)
    add_subdirectory(docs)
endif(BUILD-DOCS)