#pragma once

#include "filter.h"
#include <algorithm>
#include <vector>
#include <initializer_list>
#include "random"
#include "optional"

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
//filter.FilterTime() consists of all non zero freq return map(key, complex_t)

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree& tree, const SplittingTree::NodePtr& cone_node, const SignalInfo& info, int64_t sparsity, IndexGenerator& delta) {
    auto filter = Filter(tree, cone_node, info);
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

std::optional<FrequencyMap> SparseFFT(const Signal& x, const SignalInfo& info, int64_t expected_sparsity, SplittingTree* parent_tree,
                       SplittingTree::NodePtr& parent_node, const FrequencyMap& known_freq, IndexGenerator& delta) {
    FrequencyMap recovered_freq;
    FrequencyMap total_freq = known_freq;
    SplittingTree tree(parent_tree, parent_node);

    while (!tree.IsEmpty() && (tree.LeavesCount() + static_cast<int>(recovered_freq.size())) <= expected_sparsity) {
        NodePtr node = tree.GetLightestNode();
        if (node->level == info.Dimensions() * CalcLog(info.SignalWidth())) {
            auto filter = Filter(tree, node, info);
            complex_t recovered = 0;
            for (const auto& freq : total_freq) {
                recovered += freq.second * filter.FilterFrequency(freq.first);
            }
            complex_t filtered = 0;
            for (const auto& value : filter.FilterTime()) {
                filtered += value.second * x.ValueAtTime(-value.first);
            }
            auto value = static_cast<double>(info.SignalSize()) * filtered - recovered;
            auto key = Key(info, node->label);
            recovered_freq[key] += value;
            total_freq[key] += value;
            tree.RemoveNode(node);
        } else {
            node->AddChildren();
            if (!ZeroTest(x, total_freq, tree, node->left, info, expected_sparsity, delta)) {
                tree.RemoveNode(node->left);
            }
            if (!ZeroTest(x, total_freq, tree, node->right, info, expected_sparsity, delta)) {
                tree.RemoveNode(node->right);
            }
        }
    }
    if ((tree.LeavesCount() + static_cast<int>(recovered_freq.size())) > expected_sparsity) {
        return std::nullopt;
    }
    return recovered_freq;
}

std::optional<FrequencyMap> TryRestore(const Signal& x, const SignalInfo& info, int64_t expected_sparsity, double sparsity_step, SplittingTree* parent_tree,
                        SplittingTree::NodePtr& parent_node, const FrequencyMap& known_freq, int rank, IndexGenerator& delta) {
    FrequencyMap recovered_freq;
    FrequencyMap total_freq = known_freq;
    SplittingTree tree(parent_tree, parent_node);

    int64_t next_sparsity = std::max<int64_t>(1, llround(expected_sparsity / sparsity_step));
    while (!tree.IsEmpty() && (tree.LeavesCount() * next_sparsity + static_cast<int>(recovered_freq.size())) <= expected_sparsity) {
        NodePtr v = tree.GetLightestNode();
        std::optional<FrequencyMap> probable_freq;

        v->AddChildren();
        std::array<NodePtr, 2> children({v->left, v->right});
        for (auto& node : children) {
            if (rank == 2) {
                probable_freq = SparseFFT(x, info, next_sparsity, &tree, node, total_freq, delta);
            } else {
                probable_freq = TryRestore(x, info, next_sparsity, sparsity_step, &tree, node, total_freq, rank - 1,
                                           delta);
            }
            if (probable_freq) {
                auto new_total_freq = MapUnion(probable_freq.value(), total_freq);
                if (!ZeroTest(x, new_total_freq, tree, node, info, expected_sparsity, delta)) {
                    tree.RemoveNode(node);
                    recovered_freq = MapUnion(recovered_freq, probable_freq.value());
                    total_freq.swap(new_total_freq);
                }
            }
        }

    }
    if ((tree.LeavesCount() * next_sparsity + static_cast<int>(recovered_freq.size())) > expected_sparsity) {
        return std::nullopt;
    }
    return recovered_freq;
}

FrequencyMap RecursiveSparseFFT(const Signal& x, const SignalInfo& info, int64_t sparsity, int rank) {
    assert(info.SignalSize() > 1);
    if (sparsity == 0) {
        return {};
    }
    IndexGenerator delta(info, 61);
    NodePtr parent(nullptr);
    std::optional<FrequencyMap> res;
    if (rank == 1) {
        res = SparseFFT(x, info, sparsity, nullptr, parent, {}, delta);
    } else {
        res = TryRestore(x, info, sparsity, pow(sparsity, 1. / rank), nullptr, parent, {}, rank, delta);
    }
    assert(res);
    return res.value();
}