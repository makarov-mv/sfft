cmake_minimum_required(VERSION 3.8.2)
project(dimension_independent_SFFT)

set(CMAKE_CXX_STANDARD 17)

add_executable(test_fft src/test_fft.cpp src/filter.h src/disfft.h)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wfloat-equal -Wfloat-conversion")