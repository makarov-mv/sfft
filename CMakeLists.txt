cmake_minimum_required(VERSION 3.8.2)
project(dimension_independent_SFFT)

set(CMAKE_CXX_STANDARD 14)
#include(${CMAKE_SOURCE_DIR}/src)

add_executable(test_fft src/test_fft.cpp src/filter.h src/disfft.h)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wfloat-equal -Wconversion")