#include "catch.hpp"
#include "fftw_utility_test.h"


TEST_CASE("FFT 4") {
    const int64_t signal_size = 4;
    const int64_t sparsity = 2;
    std::vector<complex_t> out(signal_size);
    out[0] = 1;
    out[2] = 1;
    REQUIRE(RunFFTWTest(out, signal_size, sparsity));
}

TEST_CASE("FFT 32") {
    const int64_t signal_size = 32;
    const int64_t sparsity = 4;
    std::vector<complex_t> out(signal_size);
    out[0] = 1;
    out[2] = 1;
    out[9] = 93;
    out[24] = complex_t{1, -0.6};
    REQUIRE(RunFFTWTest(out, signal_size, sparsity));
}

TEST_CASE("FFT 128 arithmetic") {
    const int64_t signal_size = 128;
    std::vector<complex_t> out(signal_size);
    for (int i = 0; i < signal_size; i += 2) {
        out[i] = (i + 1.) + (3. * i - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, signal_size, 64));
}

TEST_CASE("FFT 1024 arithmetic") {
    const int64_t signal_size = 1024;
    std::vector<complex_t> out(signal_size);
    for (int i = 0; i < signal_size; i += 2) {
        out[i] = (i + 1.) + (3. * i - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, signal_size, 512));
}

TEST_CASE("FFT 1024 * 1024 sparse arithmetic") {
    const int64_t signal_size = 1024 * 1024;
    const int64_t sparsity = 16;
    std::vector<complex_t> out(signal_size);
    for (int i = 0; i < signal_size; i += signal_size / sparsity) {
        out[i] = (i + 1.) + (3. * i - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, signal_size, sparsity));
} // 1e7

TEST_CASE("FFT 1024 sparse arithmetic") {
    const int64_t signal_size = 1024;
    const int64_t sparsity = 32;
    std::vector<complex_t> out(signal_size);
    for (int i = 0; i < signal_size; i += signal_size / sparsity) {
        int id = i % 3;
        out[i] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, signal_size, sparsity));
}

TEST_CASE("FFT 1024 sparse random") {
    const int64_t signal_size = 1024;
    const int64_t sparsity = 32;
    std::vector<complex_t> out(signal_size);

    for (int i = 0; i < signal_size; i += signal_size / sparsity) {
        size_t j = random() % signal_size;
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, signal_size, sparsity));
}

std::vector<complex_t> CalcGaussian(int64_t signal_size) {
    std::vector<complex_t> res(signal_size, 0);
    for (int j = 0; j < signal_size; ++j) {
        for (int i = -30; i <= 30; ++i) {
            double pow = j / static_cast<double>(signal_size) + i;
            pow = pow * pow;
            res[j] += std::exp(-signal_size * PI * pow);
        }
    }
    return res;
}

TEST_CASE("FFT 128 Gaussian") {
    int64_t signal_size = 128;
    auto out = CalcGaussian(signal_size);
    int cnt = 0;
    for (auto& c : out) {
        cnt += NonZero(c);
    }
    REQUIRE(RunFFTWTest(out, signal_size, cnt, CheckEqualOperator(1e-5)));
}

TEST_CASE("FFT 1231123 sparse random") {
    const int64_t signal_size = 1024 * 1024;
    const int64_t sparsity = 16;
    std::vector<complex_t> out(signal_size);

    for (int i = 0; i < signal_size; i += signal_size / sparsity) {
        size_t j = random() % signal_size;
        int id = i % 9;
        out[j] = (id + 1.) + (3. * id - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, signal_size, sparsity));
}