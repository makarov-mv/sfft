#include "catch.hpp"
#include "utility_test_fftw.h"


TEST_CASE("FFT 4") {
    SignalInfo info{1, 4};
    const int64_t sparsity = 2;
    std::vector<complex_t> out(info.SignalSize());
    out[0] = 1;
    out[2] = 1;
    REQUIRE(RunFFTWTest(out, info, sparsity));
}

TEST_CASE("FFT 32") {
    SignalInfo info{1, 32};
    const int64_t sparsity = 4;
    std::vector<complex_t> out(info.SignalSize());
    out[0] = 1;
    out[2] = 1;
    out[9] = 93;
    out[24] = complex_t{1, -0.6};
    REQUIRE(RunFFTWTest(out, info, sparsity));
}

TEST_CASE("FFT 128 arithmetic") {
    SignalInfo info{1, 128};
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += 2) {
        out[i] = (i + 1.) + (3. * i - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, 64));
}

TEST_CASE("FFT 1024 arithmetic") {
    SignalInfo info{1, 1024};
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += 2) {
        out[i] = (i + 1.) + (3. * i - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, 512));
}

TEST_CASE("FFT 1024 * 1024 sparse arithmetic") {
    SignalInfo info{1, 1024 * 1024};
    const int64_t sparsity = 16;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = (i + 1.) + (3. * i - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity));
} // 1e7

TEST_CASE("FFT 1024 sparse arithmetic") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        int id = i % 3;
        out[i] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity));
}

TEST_CASE("FFT 1024 sparse random") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity));
}

TEST_CASE("FFT 2d 16") {
    SignalInfo info{2, 4};
    const int64_t sparsity = 3;
    std::vector<complex_t> out(info.SignalSize());
    out[0] = 1;
    out[2] = 1;
    out[9] = 93;
    REQUIRE(RunFFTWTest(out, info, sparsity));
}

TEST_CASE("FFT 2d 1024") {
    SignalInfo info{2, 32};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity));
}

TEST_CASE("FFT 3d 32768") {
    SignalInfo info{3, 32};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity));
}

TEST_CASE("FFT 4d 1048576 random") {
    SignalInfo info{4, 32};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        size_t j = random() % info.SignalSize();
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity));
}

TEST_CASE("FFT 4d 1048576 comb") {
    SignalInfo info{4, 32};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = (i + 1.) + (3. * i - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity));
}
