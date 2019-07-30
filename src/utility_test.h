#pragma once

#include "disfft.h"

using Node = SplittingTree::Node;
using NodePtr = SplittingTree::NodePtr;

bool CheckEqual(complex_t a, complex_t b) {
    return !NonZero(a - b);
}

std::vector<complex_t> GetSignalFromMap(const FrequencyMap& recovered_freq, const SignalInfo& info) {
    std::vector<complex_t> signal(info.SignalSize());
    for (int64_t i = 0; i < info.SignalSize(); ++i) {
        auto it = recovered_freq.find(Key{info, i});
        if (it != recovered_freq.end()) {
            signal[i] = it->second;
        } else {
            signal[i] = 0.;
        }
    }
    return signal;
}

bool RunSFFT(const DataSignal& x, const SignalInfo& info, int64_t sparsity, const std::vector<complex_t>& desired, int rank = 1, int64_t seed = 61, TransformSettings settings = {}) {
    assert(info.SignalSize() == static_cast<int64_t>(desired.size()));
    assert(sparsity <= info.SignalSize());

    FrequencyMap frequency = RecursiveSparseFFT(x, info, sparsity, rank, seed, settings);
    auto result = GetSignalFromMap(frequency, info);

    return std::equal(desired.begin(), desired.end(), result.begin(), CheckEqual);
}

bool RunSFFT(const SignalInfo& info, int64_t sparsity, const std::vector<complex_t>& data, const std::vector<complex_t>& desired, int rank = 1, int64_t seed = 61, TransformSettings settings = {}) {
    assert(info.SignalSize() == static_cast<int64_t>(data.size()));
    assert(info.SignalSize() == static_cast<int64_t>(desired.size()));
    assert(sparsity <= info.SignalSize());
    DataSignal x(info, data.data());

    return RunSFFT(x, info, sparsity, desired, rank, seed, settings);
}

using namespace std::complex_literals;
