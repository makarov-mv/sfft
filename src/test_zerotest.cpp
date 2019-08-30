#include "catch.hpp"
#include "utility_test_fftw.h"

TEST_CASE("ZeroTest 1024 32 comb") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i + 1] = CalcKernel(i * sparsity / info.SignalWidth(), sparsity);
    }
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());
    for (int rank = 1; rank <= 3; ++rank) {
        TransformSettings settings;
        std::vector<double> koefs = {1, 0.8, 0.5, 0.3, 0.1};
        for (auto koef : koefs) {
            settings.use_comb = true;
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
        out[i + 1] = CalcKernel(i * sparsity / info.SignalWidth(), sparsity);
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

TEST_CASE("ZeroTest 1024 32 rand phase") {
    SignalInfo info{1, 1024};
    const int64_t sparsity = 100;
    std::vector<complex_t> out(info.SignalSize());
    std::mt19937_64 gen(1231);
    std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);

    auto runner = FFTWRunner(info, FFTW_BACKWARD);

    for (int rank = 1; rank <= 1; ++rank) {
        TransformSettings settings;
        std::vector<int> koefs = {1, 2, 4};
        for (auto koef : koefs) {
            settings.use_comb = false;
            settings.assume_random_phase = true;
            settings.random_phase_sparsity_koef = koef;
            const int max_iter = 100;
            int ok = 0;
            for (int i = 0; i < max_iter; ++i) {
                out.assign(out.size(), 0);
                for (int j = 0; j < sparsity; ++j) {
                    out[dist(gen)] = 1;
                }
                auto in = runner.Run(out);
                auto x = DataSignal(info, in.data());
                ok += RunSFFT(x, info, sparsity, out, rank, i, settings);
            }
            WARN("rank: " << rank << ", koef: " << koef << ", " << static_cast<double>(ok) / max_iter << " successful, " << ok << "/"
                          << max_iter);
        }
    }
}

TEST_CASE("ZeroTest 3 128 32 combined for rand phase alg") {
    SignalInfo info{3, 1 << 7};
    const int64_t sparsity = 32;
    assert((sparsity & (sparsity - 1)) == 0);
    std::vector<complex_t> out(info.SignalSize());
    std::mt19937_64 gen(1231);
    std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);

    auto runner = FFTWRunner(info, FFTW_BACKWARD);

    for (int rank = 1; rank <= 1; ++rank) {
        TransformSettings settings;
        std::vector<int> koefs = {1, 2, 4, 6};
        for (auto koef : koefs) {
            settings.use_comb = true;
            settings.assume_random_phase = true;
            settings.random_phase_sparsity_koef = koef;
            const int max_iter = 100;
            int ok = 0;
            for (int i = 0; i < max_iter; ++i) {
                out.assign(out.size(), 0);
                for (int j = 0; j < sparsity / 2; ++j) {
                    out[dist(gen)] = 1;
                }
                for (int j = 0; j < info.SignalSize(); j += info.SignalSize() / (sparsity / 2)) {
                    out[j + 1] = CalcKernel(i * (sparsity / 2) / info.SignalWidth(), sparsity / 2);
                }

                auto in = runner.Run(out);
                auto x = DataSignal(info, in.data());
                ok += RunSFFT(x, info, sparsity, out, rank, i, settings);
            }
            WARN("rank: " << rank << ", koef: " << koef << ", " << static_cast<double>(ok) / max_iter << " successful, " << ok << "/"
                          << max_iter);
        }
    }
}

TEST_CASE("ZeroTest 3 128 32 combined") {
    SignalInfo info{3, 1 << 7};
    const int64_t sparsity = 32;
    std::vector<complex_t> out(info.SignalSize());
    std::mt19937_64 gen(1231);
    std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);

    auto runner = FFTWRunner(info, FFTW_BACKWARD);

    for (int rank = 1; rank <= 4; ++rank) {
        TransformSettings settings;
        std::vector<double> koefs = {0.03};
        for (auto koef : koefs) {
            settings.use_comb = false;
            settings.zero_test_koef = koef;
            const int max_iter = 100;
            int ok = 0;
            for (int i = 0; i < max_iter; ++i) {
                out.assign(out.size(), 0);
                for (int j = 0; j < sparsity; ++j) {
                    out[dist(gen)] = 1;
                }
//                for (int j = 0; j < info.SignalSize(); j += info.SignalSize() / (sparsity / 2)) {
//                    out[j] = 1;
//                }
                auto in = runner.Run(out);
                auto x = DataSignal(info, in.data());
                ok += RunSFFT(x, info, sparsity, out, rank, i, settings);
            }
            WARN("rank: " << rank << ", koef: " << koef << ", " << static_cast<double>(ok) / max_iter << " successful, " << ok << "/"
                          << max_iter);
        }
    }
}

TEST_CASE("ZeroTest 2 1024 64 4D comb randphase") {
    SignalInfo info{4, 1 << 5};
    const int64_t sparsity = 64;
    int64_t root = int(sqrt(sparsity));
    assert(root * root == sparsity);
    assert(info.SignalWidth() % root == 0);
    std::mt19937_64 gen(123423);
    std::vector<complex_t> out(info.SignalSize());
    std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);

    auto runner = FFTWRunner(info, FFTW_BACKWARD);

    for (int rank = 1; rank <= 1; ++rank) {
        TransformSettings settings;
        std::vector<int> koefs = {1, 2, 4, 6, 16};
        for (auto koef : koefs) {
            settings.use_comb = true;
            settings.assume_random_phase = true;
            settings.random_phase_sparsity_koef = koef;
            const int max_iter = 30;
            int ok = 0;
            for (int i = 0; i < max_iter; ++i) {
                out.assign(out.size(), 0);
                std::vector<complex_t> out(info.SignalSize());
                int64_t step = info.SignalWidth() * info.SignalWidth() * info.SignalWidth();
                for (int i = 0; i < info.SignalWidth(); i += info.SignalWidth() / root) {
                    for (int j = 0; j < info.SignalWidth(); j += info.SignalWidth() / root) {
                        out[i * step + j + 1] = CalcKernel((i + j) * root / info.SignalWidth(), root);
                    }
                }
                const int rand_supp = 64;
                for (int i = 0; i < rand_supp; ++i) {
                    out[dist(gen)] = 1;
                }
                auto in = runner.Run(out);
                auto x = DataSignal(info, in.data());
                ok += RunSFFT(x, info, sparsity + rand_supp, out, rank, i, settings);
            }
            WARN("rank: " << rank << ", koef: " << koef << ", " << static_cast<double>(ok) / max_iter << " successful, " << ok << "/"
                          << max_iter);
        }
    }
}

TEST_CASE("ZeroTest 2 1024 64 2D comb") {
    SignalInfo info{2, 1 << 10};
    const int64_t sparsity = 64;
    int64_t root = int(sqrt(sparsity));
    assert(root * root == sparsity);
    assert(info.SignalWidth() % root == 0);
    std::vector<complex_t> out(info.SignalSize());

    auto runner = FFTWRunner(info, FFTW_BACKWARD);

    for (int rank = 1; rank <= 4; ++rank) {
        TransformSettings settings;
        std::vector<double> koefs = {1, 0.5, 0.1};
        for (auto koef : koefs) {
            settings.use_comb = true;
            settings.zero_test_koef = koef;
            const int max_iter = 100;
            int ok = 0;
            for (int i = 0; i < max_iter; ++i) {
                out.assign(out.size(), 0);
                std::vector<complex_t> out(info.SignalSize());
                for (int i = 0; i < info.SignalWidth(); i += info.SignalWidth() / root) {
                    for (int j = 0; j < info.SignalWidth(); j += info.SignalWidth() / root) {
                        out[i * info.SignalWidth() + j + 1] = CalcKernel((i + j) * root / info.SignalWidth(), root);
                    }
                }
                auto in = runner.Run(out);
                auto x = DataSignal(info, in.data());
                ok += RunSFFT(x, info, sparsity, out, rank, i, settings);
            }
            WARN("rank: " << rank << ", koef: " << koef << ", " << static_cast<double>(ok) / max_iter << " successful, " << ok << "/"
                          << max_iter);
        }
    }
}