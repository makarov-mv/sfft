cmake_minimum_required(VERSION 3.8)
project(dimension_independent_SFFT)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wfloat-equal -Wfloat-conversion -g")

set(CMAKE_CXX_FLAGS_RELEASE "-march=native -fopenmp ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
set(CXX_FILES src/arithmetics.cpp)

add_executable(test_fft src/test_main.cpp src/test_fft.cpp ${CXX_FILES})
target_link_libraries(test_fft fftw3)
add_executable(test_fft_fftw src/test_main.cpp src/test_fft_fftw.cpp ${CXX_FILES})
target_link_libraries(test_fft_fftw fftw3)
add_executable(test_fft_fftw_rank src/test_main.cpp src/test_fft_fftw_rank.cpp ${CXX_FILES})
target_link_libraries(test_fft_fftw_rank fftw3)
add_executable(bench src/test_main.cpp src/bench.cpp ${CXX_FILES})
target_link_libraries(bench fftw3)
add_executable(test_zerotest src/test_main.cpp src/test_zerotest.cpp ${CXX_FILES})
target_link_libraries(test_zerotest fftw3)
add_executable(test_projrec src/test_main.cpp src/test_projrec.cpp ${CXX_FILES})
target_link_libraries(test_projrec fftw3)
add_executable(measure_run src/measure_run.cpp ${CXX_FILES})
target_link_libraries(measure_run fftw3)
