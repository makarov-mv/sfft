#pragma once

#include "complex"
#include "assert.h"
#include "memory"
#include "unordered_map"
#include "vector"
#include "random"

using complex_t  = std::complex<double>;

struct TreeLabel {
    int length{0};
    uint64_t label{0};

    TreeLabel CalcParent() const {
        assert(length > 0);
        return {length - 1, label & ((1ull << static_cast<uint64_t>(length - 1)) - 1)};
    }

    TreeLabel CalcLeft() const {
        assert(length > 0);
        return {length + 1, label + (1ull << static_cast<uint64_t>(length + 1))};
    }

    TreeLabel CalcRight() const {
        assert(length > 0);
        return {length + 1, label};
    }
};

class SplittingTree {
public:
    struct Node;
    using NodePtr = std::shared_ptr<Node>;
    struct Node {
        TreeLabel label;
        NodePtr left{nullptr};
        NodePtr right{nullptr};
        std::weak_ptr<Node> parent;
    };

    SplittingTree(): root_(std::make_shared<Node>()) {
    }

private:

    std::shared_ptr<Node> root_;
};

using FrequencyMap = std::unordered_map<uint64_t, complex_t>;

class Signal {
public:
    virtual complex_t ValueAtTime(int64_t) const;
};

class SignalVector: public Signal {
public:
    SignalVector(std::vector<complex_t> v): values_(v) {
    }

    complex_t ValueAtTime(int64_t t) const override {
        return values_[t];
    }

private:
    std::vector<complex_t> values_;
};

class IndexGenerator {
public:
    IndexGenerator(uint64_t signal_size, uint64_t seed): rand_gen_(seed), index_gen_(0, signal_size) {
    }

    uint64_t Next() {
        return index_gen_(rand_gen_);
    }

private:
    std::mt19937_64 rand_gen_;
    std::uniform_int_distribution<uint64_t> index_gen_;
    uint64_t
};

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree& tree, const SplittingTree::NodePtr& cone_node, uint64_t signal_size, IndexGenerator& delta) {

}