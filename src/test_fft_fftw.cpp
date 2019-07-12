#include "catch.hpp"
#include "utility_test.h"
#include "fftw3.h"

std::vector<fftw_complex> MakeFFTWVector(std::vector<complex_t> x) {
    std::vector<fftw_complex> res(x.size());
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        res[i][0] = x[i].real();
        res[i][1] = x[i].imag();
    }
}

TEST_CASE("SFFT 16") {
    const int64_t signal_size = 16;
    std::vector<fftw_complex> in(signal_size), out(signal_size);
    fftw_plan plan;
    plan = fftw_plan_dft_1d(signal_size, in.data(), out.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
}
