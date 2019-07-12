#pragma once

#include "filter.h"
#include <algorithm>
#include "random"

using FrequencyMap = std::unordered_map<int64_t, complex_t>;
using NodePtr = SplittingTree::NodePtr;
using Node = SplittingTree::Node;

class Signal {
public:
    virtual ~Signal() = default;

    virtual complex_t ValueAtTime(int64_t) const = 0;
};

class DataSignal: public Signal {
public:
    DataSignal(int64_t signal_size, const complex_t* v): signal_size_(signal_size), values_(v) {
    }

    complex_t ValueAtTime(int64_t t) const override {
        t %= signal_size_;
        if (t < 0) {
            t += signal_size_;
        }
        return values_[t];
    }

private:
    const int64_t signal_size_;
    const complex_t* values_;
};

class IndexGenerator {
public:
    IndexGenerator(int64_t signal_size, int64_t seed): rand_gen_(seed), index_gen_(0, signal_size) {
    }

    int64_t Next() {
        return index_gen_(rand_gen_);
    }

private:
    std::mt19937_64 rand_gen_;
    std::uniform_int_distribution<int64_t> index_gen_;
};

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree::NodePtr& cone_node, int64_t signal_size, int64_t sparsity, IndexGenerator& delta) {
    auto filter = Filter(cone_node, signal_size);
    int64_t max_iters = std::max<int64_t>(llround(2 * sparsity * log2(sparsity) * log2(sparsity) * log2(signal_size)), 2);
    for (int64_t i = 0; i < max_iters; ++i) {
        auto time = delta.Next();
        complex_t recovered_at_time = 0;
        for (auto freq: recovered_freq) {
            recovered_at_time += CalcKernel(freq.first * time, signal_size) * freq.second * filter.FilterFrequency(freq.first) / complex_t(signal_size, 0);
        }
        complex_t filtered_at_time = 0;
        for (auto value: filter.FilterTime()) {
            filtered_at_time += value.second * x.ValueAtTime(time - value.first);
        }
        if (NonZero(filtered_at_time - recovered_at_time)) {
            return true;
        }
    }
    return false;
}

FrequencyMap SparseFFT(const Signal& x, int64_t signal_size, int64_t sparsity) {
    assert(signal_size > 1);
    SplittingTree tree{};
    FrequencyMap recovered_freq;
    IndexGenerator delta{signal_size, 61};

    while (tree.IsNonEmpty()) {
        NodePtr node = tree.GetLightestNode();
        if (node->level == CalcLog(signal_size)) {
            auto filter = Filter(node, signal_size);
            complex_t recovered = 0;
            for (auto freq : recovered_freq) {
                recovered += freq.second * filter.FilterFrequency(freq.first);
            }
            complex_t filtered = 0;
            for (auto value : filter.FilterTime()) {
                filtered += value.second * x.ValueAtTime((signal_size - value.first) % signal_size);
            }
            recovered_freq[node->label] += signal_size * 1. * filtered - recovered;
            tree.RemoveNode(node);
        } else {
            node->AddChildren();
            if (!ZeroTest(x, recovered_freq, node->left, signal_size, sparsity, delta)) {
                node->left = nullptr;
            }
            if (!ZeroTest(x, recovered_freq, node->right, signal_size, sparsity, delta)) {
                node->right = nullptr;
            }
        }
    }
    return recovered_freq;
}