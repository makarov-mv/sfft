cmake_minimum_required(VERSION 3.8.2)
project(sFFT)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wfloat-equal")
set(SRC computefourier.cc fftw.cc filters.cc parameters.cc timer.cc utils.cc plot.cc)

add_executable(experiment experiment.cc ${SRC})
target_link_libraries(experiment fftw3)
add_executable(generate_graphs generate_graphs.cc ${SRC})
target_link_libraries(generate_graphs fftw3)

