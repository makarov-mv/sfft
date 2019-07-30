#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include "utility_test_fftw.h"

TEST_CASE("benchmark 1024 comb 32") {
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

TEST_CASE("benchmark 1024 comb 32 no preempt") {
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
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 1, i, false); });
        };
    BENCHMARK_ADVANCED("rank 2")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 2, i, false); });
        };
    BENCHMARK_ADVANCED("rank 3")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 3, i, false); });
        };
    BENCHMARK_ADVANCED("rank 4")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 4, i, false); });
        };
}

TEST_CASE("benchmark 1024 random 32") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
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

TEST_CASE("benchmark 1024 random 32 no preempt") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());

    BENCHMARK_ADVANCED("rank 1")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 1, i, false); });
        };
    BENCHMARK_ADVANCED("rank 2")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 2, i, false); });
        };
    BENCHMARK_ADVANCED("rank 3")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 3, i, false); });
        };
    BENCHMARK_ADVANCED("rank 4")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 4, i, false); });
        };
}

TEST_CASE("benchmark 1024 random 32 no preempt") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());

    BENCHMARK_ADVANCED("rank 1")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 1, i, false); });
        };
    BENCHMARK_ADVANCED("rank 2")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 2, i, false); });
        };
    BENCHMARK_ADVANCED("rank 3")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 3, i, false); });
        };
    BENCHMARK_ADVANCED("rank 4")(Catch::Benchmark::Chronometer meter) {
            meter.measure([&](int i) { return RecursiveSparseFFT(x, info, sparsity, 4, i, false); });
        };
}

TEST_CASE("benchmark 4d 1048576 random 32") {
    SignalInfo info{4, 32};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
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