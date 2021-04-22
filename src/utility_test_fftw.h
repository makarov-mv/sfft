#pragma once

#include "utility_test.h"
#include "fftwrunner.h"

bool RunFFTWTest(const std::vector<complex_t>& out, const SignalInfo& info, int64_t sparsity, int rank = 1, int64_t seed = 61, TransformSettings settings = {}) {
    assert(info.SignalSize() == static_cast<int64_t>(out.size()));
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());
    auto result = GetSignalFromMap(RecursiveSparseFFT(x, info, sparsity, rank, seed, settings), info);
    return std::equal(out.begin(), out.end(), result.begin(), CheckEqual);
}