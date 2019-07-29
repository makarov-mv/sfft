#include "catch.hpp"
#include "utility_test.h"
#include "fftw3.h"

class FFTWRunner {
public:
    FFTWRunner(const SignalInfo& info, int direction): info_(info), ranks_(info.Dimensions(), info_.SignalWidth()), direction_(direction) {
        in_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * info.SignalSize());
        out_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * info.SignalSize());
        plan_ = fftw_plan_dft(info.Dimensions(), ranks_.data(), in_, out_, direction, FFTW_ESTIMATE);
    }

    ~FFTWRunner() {
        fftw_destroy_plan(plan_);
        fftw_free(in_);
        fftw_free(out_);
        fftw_cleanup();
    }

    std::vector<complex_t> Run(const std::vector<complex_t>& x) {
        assert(static_cast<int64_t>(x.size()) == info_.SignalSize());
        for (int i = 0; i < info_.SignalSize(); ++i) {
            in_[i][0] = x[i].real();
            in_[i][1] = x[i].imag();
        }

        fftw_execute(plan_);

        auto res = std::vector<complex_t>(info_.SignalSize());
        for (int i = 0; i < info_.SignalSize(); ++i) {
            res[i] = {out_[i][0], out_[i][1]};
            if (direction_ == FFTW_BACKWARD) {
                res[i] /= info_.SignalSize();
            }
        }
        return res;
    }

private:
    SignalInfo info_;
    const std::vector<int> ranks_;
    int direction_;
    fftw_complex* in_;
    fftw_complex* out_;
    fftw_plan plan_;
};

bool RunFFTWTest(const std::vector<complex_t>& out, const SignalInfo& info, int64_t sparsity) {
    assert(info.SignalSize() == static_cast<int64_t>(out.size()));
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());
    auto result = GetSignalFromMap(SparseFFT(x, info, sparsity), info);
    return std::equal(out.begin(), out.end(), result.begin(), CheckEqual);
}


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

TEST_CASE("FFT 4d 1048576 comb") {
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

TEST_CASE("FFT 4d 1048576 random") {
    SignalInfo info{4, 32};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = (i + 1.) + (3. * i - 1.) * 1.i;
    }
    REQUIRE(RunFFTWTest(out, info, sparsity));
}