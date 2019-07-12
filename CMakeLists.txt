cmake_minimum_required(VERSION 3.8.2)
project(dimension_independent_SFFT)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wfloat-equal -Wfloat-conversion")

add_executable(test_fft src/test_main.cpp src/test_fft.cpp)
add_executable(test_fft_fftw src/test_main.cpp src/test_fft_fftw.cpp)
target_link_libraries(test_fft_fftw fftw3)
