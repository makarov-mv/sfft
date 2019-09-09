#pragma once

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
//        fftw_cleanup();
    }

    std::vector<complex_t> Run(const complex_t* x) {
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

    std::vector<complex_t> Run(const std::vector<complex_t>& x) {
        assert(static_cast<int64_t>(x.size()) == info_.SignalSize());
        return Run(x.data());
    }

private:
    SignalInfo info_;
    const std::vector<int> ranks_;
    int direction_;
    fftw_complex* in_;
    fftw_complex* out_;
    fftw_plan plan_;
};

bool RunFFTWTest(const std::vector<complex_t>& out, const SignalInfo& info, int64_t sparsity, int rank = 1, int64_t seed = 61, TransformSettings settings = {}) {
    assert(info.SignalSize() == static_cast<int64_t>(out.size()));
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());
    auto result = GetSignalFromMap(RecursiveSparseFFT(x, info, sparsity, rank, seed, settings), info);
    return std::equal(out.begin(), out.end(), result.begin(), CheckEqual);
}