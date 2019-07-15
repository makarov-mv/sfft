#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include "catch.hpp"
#include "fftw_utility_test.h"

//TEST_CASE("sfft 1024 32") {
//    const int64_t signal_size = 1024;
//    const int64_t sparsity = 32;
//    std::vector<complex_t> out(signal_size);
//    auto runner = FFTWRunner(signal_size, FFTW_BACKWARD);
//    {
//        for (int i = 0; i < signal_size; i += signal_size / sparsity) {
//            int id = i % 3;
//            out[i] = (id + 1.) + (3. * id - 1.) * 1.i;
//        }
//
//        auto in = runner.Run(out);
//        auto x = DataSignal(signal_size, in.data());
//
//        BENCHMARK("benchmark 1") {
//            return SparseFFT(x, signal_size, sparsity);
//        };
//    }
//    {
//        out.assign(signal_size, 0);
//        for (int i = 0; i < signal_size; i += signal_size / sparsity) {
//            out[i] = 1.;
//        }
//
//        auto in = runner.Run(out);
//        auto x = DataSignal(signal_size, in.data());
//
//        BENCHMARK("benchmark 2") {
//            return SparseFFT(x, signal_size, sparsity);
//        };
//    }
//}

complex_t f1(double power, double base) {
    double phi = 2 * PI * power / base;
    return {cos(phi), sin(phi)};
}

TEST_CASE("kernel") {
    BENCHMARK("f1") {
        return f1(10, 128);
    };
    
}