#pragma once

#include "filter.h"
#include <algorithm>
#include <vector>
#include <initializer_list>
#include "random"

class Signal {
public:
    virtual ~Signal() = default;

    virtual complex_t ValueAtTime(const Key& key) const = 0;
};

class DataSignal: public Signal {
public:
    DataSignal(const SignalInfo& info, const complex_t* v):
            info_(info),
            values_(v) {
    }

    complex_t ValueAtTime(const Key& key) const override {
        assert(key.GetSignalInfo() == info_);
        return values_[key.Flatten()];
    }

private:
    const SignalInfo info_;
    const complex_t* values_;
};

class IndexGenerator {
public:
    IndexGenerator(const SignalInfo& info, int64_t seed):
        info_(info), rand_gen_(seed), index_gen_(0, info.SignalSize()) {
    }

    Key Next() {
        return Key(info_, index_gen_(rand_gen_));
    }

private:
    SignalInfo info_;
    std::mt19937_64 rand_gen_;
    std::uniform_int_distribution<int64_t> index_gen_;
};

//filter.FilterFrequency(const Key& key)
//filter.FilterTime() consists of all no zero freq return map(key, complex_t)

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree::NodePtr& cone_node, const SignalInfo& info, int64_t sparsity, IndexGenerator& delta) {
    auto filter = Filter(cone_node, info);
    int64_t max_iters = std::max<int64_t>(llround(2 * sparsity * log2(info.SignalSize())), 2); // check
    for (int64_t i = 0; i < max_iters; ++i) {
        auto time = delta.Next();
        complex_t recovered_at_time = 0;
        for (const auto& freq: recovered_freq) {
            recovered_at_time += CalcKernel(freq.first * time, info.SignalWidth()) * freq.second * filter.FilterFrequency(freq.first);
        }
        recovered_at_time /= static_cast<double>(info.SignalSize());
        complex_t filtered_at_time = 0;
        for (const auto& value: filter.FilterTime()) {
            filtered_at_time += value.second * x.ValueAtTime(time - value.first);
        }
        if (NonZero(filtered_at_time - recovered_at_time)) {
            return true;
        }
    }
    return false;
}

FrequencyMap SparseFFT(const Signal& x, const SignalInfo& info, int64_t sparsity) {
    assert(info.SignalSize() > 1);
    SplittingTree tree{};
    FrequencyMap recovered_freq;
    IndexGenerator delta{info, 61};

    while (!tree.IsEmpty()) {
        NodePtr node = tree.GetLightestNode();
        if (node->level == info.Dimensions() * CalcLog(info.SignalWidth())) {
            auto filter = Filter(node, info);
            complex_t recovered = 0;
            for (const auto& freq : recovered_freq) {
                recovered += freq.second * filter.FilterFrequency(freq.first);
            }
            complex_t filtered = 0;
            for (const auto& value : filter.FilterTime()) {
                filtered += value.second * x.ValueAtTime(-value.first);
            }
            recovered_freq[Key(info, node->label)] += static_cast<double>(info.SignalSize()) * filtered - recovered;
            tree.RemoveNode(node);
        } else {
            node->AddChildren();
            if (!ZeroTest(x, recovered_freq, node->left, info, sparsity, delta)) {
                tree.RemoveNode(node->left);
            }
            if (!ZeroTest(x, recovered_freq, node->right, info, sparsity, delta)) {
                tree.RemoveNode(node->right);
            }
        }
    }
    return recovered_freq;
}