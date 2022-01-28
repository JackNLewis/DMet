cmake_minimum_required(VERSION 3.20)
project(DMet)

set(CMAKE_CXX_STANDARD 14)

include_directories(include/)
add_library(DMet SHARED src/DMet/DistribDistances.cpp src/DMet/PointDistances.cpp src/DMet/utils.cpp)


##include gmp headers and mpfr
include_directories("/usr/local/include/")

add_subdirectory(test)