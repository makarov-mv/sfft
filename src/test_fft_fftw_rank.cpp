#include "catch.hpp"
#include "utility_test_fftw.h"

TEST_CASE("FFT 1024 rank 2 comb") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = 1;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity, 2));
}

TEST_CASE("FFT 1024 rank 3 comb") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = 1;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity, 3));
}

TEST_CASE("FFT 1024 rank 2 random") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity, 2));
}

TEST_CASE("FFT 1024 rank 3 random") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity, 3));
}

TEST_CASE("FFT 3d rank 3 32768") {
    SignalInfo info{3, 32};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity, 3));
}

TEST_CASE("FFT 4d rank 4 1048576 random") {
    SignalInfo info{4, 32};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity, 4));
}