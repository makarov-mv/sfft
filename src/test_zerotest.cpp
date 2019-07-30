#include "catch.hpp"
#include "utility_test_fftw.h"

TEST_CASE("ZeroTest 1024 32 comb") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = 1;
    }
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());
    for (int rank = 1; rank <= 3; ++rank) {
        TransformSettings settings;
        std::vector<double> koefs = {1, 0.8, 0.6, 0.4, 0.2};
        for (auto koef : koefs) {
            settings.zero_test_koef = koef;
            const int max_iter = 100;
            int ok = 0;
            for (int i = 0; i < max_iter; ++i) {
                ok += RunSFFT(x, info, sparsity, out, rank, i, settings);
            }
            WARN("rank: " << rank << ", koef: " << koef << ", " << static_cast<double>(ok) / max_iter << " successful, " << ok << "/"
                          << max_iter);
        }
    }
}

TEST_CASE("ZeroTest 1024 * 1024 32 comb") {
    SignalInfo info{1, 1024 * 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = 1;
    }
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());

    TransformSettings settings;
    std::vector<double> koefs = {1, 0.8, 0.6, 0.4, 0.2};
    for (int rank = 1; rank <= 3; ++rank) {
        for (auto koef : koefs) {
            settings.zero_test_koef = koef;
            const int max_iter = 100;
            int ok = 0;
            for (int i = 0; i < max_iter; ++i) {
                ok += RunSFFT(x, info, sparsity, out, 1, i, settings);
            }
            WARN("rank: " << rank << ", koef: " << koef << ", " << static_cast<double>(ok) / max_iter << " successful, " << ok << "/"
                          << max_iter);
        }
    }
}