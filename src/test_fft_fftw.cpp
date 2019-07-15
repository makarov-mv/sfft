#include "catch.hpp"
#include "utility_test.h"
#include "fftw3.h"

class FFTWRunner {
public:
    FFTWRunner(int64_t signal_size, int direction): signal_size_(signal_size), direction_(direction) {
        in_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * signal_size);
        out_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * signal_size);
        plan_ = fftw_plan_dft_1d(signal_size, in_, out_, direction, FFTW_ESTIMATE);
    }

    ~FFTWRunner() {
        fftw_destroy_plan(plan_);
        fftw_free(in_);
        fftw_free(out_);
        fftw_cleanup();
    }

    std::vector<complex_t> Run(const std::vector<complex_t>& x) {
        assert(static_cast<int64_t>(x.size()) == signal_size_);
        for (int i = 0; i < signal_size_; ++i) {
            in_[i][0] = x[i].real();
            in_[i][1] = x[i].imag();
        }

        fftw_execute(plan_);

        auto res = std::vector<complex_t>(signal_size_);
        for (int i = 0; i < signal_size_; ++i) {
            res[i] = {out_[i][0], out_[i][1]};
            if (direction_ == FFTW_BACKWARD) {
                res[i] /= signal_size_;
            }
        }
        return res;
    }

private:
    int64_t signal_size_;
    int direction_;
    fftw_complex* in_;
    fftw_complex* out_;
    fftw_plan plan_;
};

bool RunFFTWTest(const std::vector<complex_t>& out, int64_t signal_size, int64_t sparsity) {
    assert(signal_size = static_cast<int64_t>(out.size()));
    auto runner = FFTWRunner(signal_size, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(signal_size, in.data());
    auto result = GetSignalFromMap(SparseFFT(x, signal_size, sparsity), signal_size);
    return std::equal(out.begin(), out.end(), result.begin(), CheckEqual);
}


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

TEST_CASE("FFT stress 128 16") {
    const int64_t signal_size = 128;
    std::vector<complex_t> out(signal_size);
    out[0] = 1;
    out[2] = 1;
    out[9] = 93;
    out[24] = complex_t{1, -0.6};
    REQUIRE(RunFFTWTest(out, signal_size, 16));
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