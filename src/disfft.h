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

struct TransformSettings {
    bool use_preemptive_tests{true};
    double zero_test_koef{1};
};

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree& tree,
              const SplittingTree::NodePtr& cone_node, const SignalInfo& info, int64_t sparsity, IndexGenerator& delta,
              const TransformSettings& settings) {
    auto filter = Filter(tree, cone_node, info);
    int64_t max_iters = std::max<int64_t>(llround(settings.zero_test_koef * sparsity * log2(info.SignalSize())), 1);
    std::vector<std::pair<Key, complex_t>> freq_precalc;
    freq_precalc.reserve(recovered_freq.size());
    for (const auto& freq: recovered_freq) {
        freq_precalc.emplace_back(freq.first, freq.second * filter.FilterFrequency(freq.first));
    }
    for (int64_t i = 0; i < max_iters; ++i) {
        auto time = delta.Next();
        complex_t recovered_at_time = 0;
        for (const auto& freq: freq_precalc) {
            recovered_at_time += CalcKernel(freq.first * time, info.SignalWidth()) * freq.second;
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
                       SplittingTree::NodePtr& parent_node, const FrequencyMap& known_freq, IndexGenerator& delta, const TransformSettings& settings) {
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
            if (!ZeroTest(x, total_freq, tree, node->left, info, expected_sparsity, delta, settings)) {
                tree.RemoveNode(node->left);
            }
            if (!ZeroTest(x, total_freq, tree, node->right, info, expected_sparsity, delta, settings)) {
                tree.RemoveNode(node->right);
            }
        }
    }
    if ((tree.LeavesCount() + static_cast<int>(recovered_freq.size())) > expected_sparsity) {
        return std::nullopt;
    }
    return recovered_freq;
}

class Restorer {
public:
    Restorer(const SplittingTree* parent_tree, SplittingTree::NodePtr& parent_node, const FrequencyMap& known_freq)
        : recovered_freq_(), total_freq_(known_freq), tree_(parent_tree, parent_node) {}

    std::optional<FrequencyMap> TryRestore(const Signal& x, const SignalInfo& info, int64_t expected_sparsity, double sparsity_step,
                                           int rank, IndexGenerator& delta, const TransformSettings& settings) {

        next_sparsity_ = std::max<int64_t>(1, llround(expected_sparsity / sparsity_step));
        auto root = tree_.GetRoot();
        SparsityTest(x, info, expected_sparsity, sparsity_step, root, rank, delta, settings);
        while (!tree_.IsEmpty() && (tree_.LeavesCount() * next_sparsity_ + static_cast<int>(recovered_freq_.size())) <= expected_sparsity) {
            NodePtr v = tree_.GetLightestNode();
            v->AddChildren();
            std::array<NodePtr, 2> children({v->left, v->right});
            bool skip_restore = false;
            if (settings.use_preemptive_tests) {
                for (auto& node : children) {
                    if (!ZeroTest(x, total_freq_, tree_, node, info, expected_sparsity, delta, settings)) {
                        skip_restore = true;
                        tree_.RemoveNode(node);
                    }
                }
            }
            if (!skip_restore) {
                for (auto& node : children) {
                    SparsityTest(x, info, expected_sparsity, sparsity_step, node, rank, delta, settings);
                }
            }
        }
        if ((tree_.LeavesCount() * next_sparsity_ + static_cast<int>(recovered_freq_.size())) > expected_sparsity) {
            return std::nullopt;
        }
        return std::move(recovered_freq_);
    }

private:
    void SparsityTest(const Signal& x, const SignalInfo& info, int64_t expected_sparsity, double sparsity_step,
                                           SplittingTree::NodePtr& node, int rank, IndexGenerator& delta, const TransformSettings& settings) {
        std::optional<FrequencyMap> probable_freq(std::nullopt);
        if (rank == 2) {
            probable_freq = SparseFFT(x, info, next_sparsity_, &tree_, node, total_freq_, delta, settings);
        } else {
            Restorer restorer(&tree_, node, total_freq_);
            probable_freq = restorer.TryRestore(x, info, next_sparsity_, sparsity_step, rank - 1,
                                                delta, settings);
        }
        if (probable_freq) {
            auto new_total_freq = MapUnion(probable_freq.value(), total_freq_);
            if (!ZeroTest(x, new_total_freq, tree_, node, info, expected_sparsity, delta, settings)) {
                tree_.RemoveNode(node);
                recovered_freq_ = MapUnion(recovered_freq_, probable_freq.value());
                total_freq_.swap(new_total_freq);
            }
        }
    }

    FrequencyMap recovered_freq_;
    FrequencyMap total_freq_;
    int64_t next_sparsity_;
    SplittingTree tree_;
};

FrequencyMap RecursiveSparseFFT(const Signal& x, const SignalInfo& info, int64_t sparsity, int rank, int64_t seed = 61, TransformSettings settings = {}) {
    assert(info.SignalSize() > 1);
    if (sparsity == 0) {
        return {};
    }
    IndexGenerator delta(info, seed);
    NodePtr parent(nullptr);
    std::optional<FrequencyMap> res;
    if (rank == 1) {
        res = SparseFFT(x, info, sparsity, nullptr, parent, {}, delta, settings);
    } else {
        Restorer restorer(nullptr, parent, {});
        res = restorer.TryRestore(x, info, sparsity, pow(sparsity, 1. / rank), rank, delta, settings);
    }
    if (!res) {
        return {};
    }
    return res.value();
}