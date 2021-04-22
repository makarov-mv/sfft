#pragma once
#include "signal.h"
#include "filter.h"
#include "fftwrunner.h"
#include "complex"

FrequencyMap RecoverByProjecting(const Signal& x, const SignalInfo& info, int projection_dim, int measures_num, int64_t transform_size, FFTWRunner& runner) {
    int64_t shift = info.LogSignalWidth() * projection_dim;
    int64_t mask = (1 << shift) - 1;

    std::vector<std::vector<complex_t>> transforms;
    Key true_i(info);
    FrequencyMap result;
    for (int t = 0; t < measures_num; ++t) {
        for (int64_t i = 0; i < transform_size; ++i) {
            true_i.SetFromFlatten(((i >> shift) << (shift + info.LogSignalWidth())) + (t << shift) + (i & mask));
            runner.SetInput(i, x.ValueAtTime(true_i));
        }
        transforms.emplace_back(std::move(runner.RunOnInput()));
    }
    for (int64_t j = 0; j < transform_size; ++j) {
        if (!NonZero(transforms[0][j])) {
            continue;
        }
        complex_t b = transforms[1][j] / transforms[0][j];
        double phase = std::arg(b);
        if (phase < 0) {
            phase += 2 * PI;
        }
        int64_t i = (std::llround(phase * info.SignalWidth() / (2 * PI))) & (info.SignalWidth() - 1);
        complex_t s = transforms[0][j];
        bool correctly_recovered = true;
        for (int t = 0; t < measures_num; ++t) {
            auto expected_val = s * CalcKernel(t * i, info.SignalWidth());
            if (NonZero(transforms[t][j] - expected_val)) {
                correctly_recovered = false;
                break;
            }
        }
        if (correctly_recovered) {
            true_i.SetFromFlatten(((j >> shift) << (shift + info.LogSignalWidth())) + (i << shift) + (j & mask));
            result.emplace(true_i, s * static_cast<double>(info.SignalWidth())); // need to rescale fsr
        }
    }
    return result;
}

FrequencyMap ProjectionFT(const Signal& x, const SignalInfo& info, int extra_measures = 6) {
    if (info.Dimensions() == 1) {
        return {}; // this case is not covered or interesting
    }
    int measures_num = extra_measures + 2;
    if (info.SignalWidth() < measures_num) {
        return {}; // don't deal with small inputs
    }
    int64_t transform_size = info.SignalSize() / info.SignalWidth();
    FFTWRunner runner(transform_size, info.Dimensions() - 1, info.SignalWidth(), FFTW_FORWARD);
    FrequencyMap result;
    for (int i = 0; i < info.Dimensions(); ++i) {
        auto call_res = RecoverByProjecting(x, info, i, measures_num, transform_size, runner);
        result.merge(call_res);
    }
    return result;
}
