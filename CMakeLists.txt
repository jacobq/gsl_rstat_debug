cmake_minimum_required(VERSION 3.10)
project(rstat_debug C)

set(CMAKE_C_STANDARD 11)

add_executable(rstat_debug main.c)

find_package(GSL REQUIRED)
target_link_libraries(rstat_debug GSL::gsl GSL::gslcblas)