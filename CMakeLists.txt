cmake_minimum_required(VERSION 3.20)
project(DMet)

set(CMAKE_CXX_STANDARD 14)

include_directories(include/)
add_library(DMet SHARED src/DMet/DistribDistances.cpp src/DMet/PointDistances.cpp src/DMet/utils.cpp)


##include gmp headers and mpfr
include_directories("/usr/local/include/")
target_link_libraries(DMet PUBLIC "/usr/local/Cellar/gmp/6.2.1_1/lib/libgmp.a")
target_link_libraries(DMet PUBLIC "/usr/local/Cellar/mpfr/4.1.0/lib/libmpfr.a")


add_subdirectory(test)