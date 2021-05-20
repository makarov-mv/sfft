#pragma once

#include "fftw3.h"
#include "assert.h"

class FFTWRunner {
public:
    FFTWRunner(int64_t signal_size, const std::vector<int>& signal_width, int direction)
    : signal_size_(signal_size), ranks_(signal_width.rbegin(), signal_width.rend()), direction_(direction) {
        in_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * signal_size);
        out_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * signal_size);
        plan_ = fftw_plan_dft(ranks_.size(), ranks_.data(), in_, out_, direction, FFTW_ESTIMATE);
    }

    FFTWRunner(const SignalInfo& info, int direction): FFTWRunner(info.SignalSize(), info.GetAllDimensions(), direction) {
    }

    ~FFTWRunner() {
        fftw_destroy_plan(plan_);
        fftw_free(in_);
        fftw_free(out_);
//        fftw_cleanup();
    }

    void SetInput(int i, complex_t value) {
        in_[i][0] = value.real();
        in_[i][1] = value.imag();
    }

    std::vector<complex_t> RunOnInput() {
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

    std::vector<complex_t> Run(const complex_t* x) {
        for (int i = 0; i < signal_size_; ++i) {
            in_[i][0] = x[i].real();
            in_[i][1] = x[i].imag();
        }

        return RunOnInput();
    }

    std::vector<complex_t> Run(const std::vector<complex_t>& x) {
        assert(static_cast<int64_t>(x.size()) == signal_size_);
        return Run(x.data());
    }

private:
    int64_t signal_size_;
    const std::vector<int> ranks_;
    int direction_;
    fftw_complex* in_;
    fftw_complex* out_;
    fftw_plan plan_;
};