#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include "utility_test_fftw.h"

TEST_CASE("benchmark 1024 comb 32") {
    // 1024 rank 2 comb
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = 1;
    }
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());

    BENCHMARK_ADVANCED("rank 1")(Catch::Benchmark::Chronometer meter) {
        meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 1, i); });
    };
    BENCHMARK_ADVANCED("rank 2")(Catch::Benchmark::Chronometer meter) {
        meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 2, i); });
    };
    BENCHMARK_ADVANCED("rank 3")(Catch::Benchmark::Chronometer meter) {
        meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 3, i); });
    };
    BENCHMARK_ADVANCED("rank 4")(Catch::Benchmark::Chronometer meter) {
        meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 4, i); });
    };
}