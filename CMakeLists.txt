cmake_minimum_required(VERSION 3.12)
project(dimension_independent_SFFT)

set(CMAKE_CXX_STANDARD 14)

add_executable(test_fft test_fft.cpp)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wfloat-equal -Wconversion")