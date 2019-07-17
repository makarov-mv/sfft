#pragma once

#include "disfft.h"

using Node = SplittingTree::Node;
using NodePtr = SplittingTree::NodePtr;

bool CheckEqual(complex_t a, complex_t b) {
    return !NonZero(a - b);
}

auto CheckEqualOperator(double eps) {
    return [eps](complex_t a, complex_t b) {
        return abs(a - b) < eps;
    };
}

std::vector<complex_t> GetSignalFromMap(const FrequencyMap& recovered_freq, int64_t signal_size) {
    std::vector<complex_t> signal(signal_size);
    for (int64_t i = 0; i < signal_size; ++i) {
        auto it = recovered_freq.find(i);
        if (it != recovered_freq.end()) {
            signal[i] = recovered_freq.at(i);
        } else {
            signal[i] = 0.;
        }
    }
    return signal;
}

bool RunSFFT(int64_t signal_size, int64_t sparsity, const std::vector<complex_t>& data, const std::vector<complex_t>& desired) {
    assert(signal_size == static_cast<int64_t>(data.size()));
    assert(signal_size == static_cast<int64_t>(desired.size()));
    assert(sparsity <= signal_size);
    DataSignal x(signal_size, data.data());

    FrequencyMap frequency = SparseFFT(x, signal_size, sparsity);
    auto result = GetSignalFromMap(frequency, signal_size);

    return std::equal(desired.begin(), desired.end(), result.begin(), CheckEqual);
}

using namespace std::complex_literals;