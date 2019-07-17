#include "utility_test.h"
#include "fftw3.h"
#include "catch.hpp"

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

template <class Equal>
bool TestSFFT(const std::vector<complex_t>& in, const std::vector<complex_t>& out, int64_t signal_size, int64_t sparsity, Equal eq) {
    assert(signal_size = static_cast<int64_t>(out.size()));
    assert(signal_size = static_cast<int64_t>(in.size()));
    auto x = DataSignal(signal_size, in.data());
    auto map = SparseFFT(x, signal_size, sparsity);
    printf("%lu\n", map.size());
    auto result = GetSignalFromMap(map, signal_size);
    return std::equal(out.begin(), out.end(), result.begin(), eq);
}

template <class Equal>
bool RunFFTWTest(const std::vector<complex_t>& out, int64_t signal_size, int64_t sparsity, Equal eq) {
    assert(signal_size = static_cast<int64_t>(out.size()));
    auto runner = FFTWRunner(signal_size, FFTW_BACKWARD);
    auto in = runner.Run(out);
    return TestSFFT(in, out, signal_size, sparsity, eq);
}

bool TestSFFT(const std::vector<complex_t>& in, const std::vector<complex_t>& out, int64_t signal_size, int64_t sparsity) {
    return TestSFFT(in, out, signal_size, sparsity, CheckEqual);
}

bool RunFFTWTest(const std::vector<complex_t>& out, int64_t signal_size, int64_t sparsity) {
    return RunFFTWTest(out, signal_size, sparsity, CheckEqual);
}