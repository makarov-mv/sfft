#pragma once

#include "complex"
#include "assert.h"
#include "memory"

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
    SplittingTree(): root_(std::make_shared<Node>()) {
    }

private:
    struct Node;
    using NodePtr = std::shared_ptr<Node>;
    struct Node {
        TreeLabel label;
        NodePtr left{nullptr};
        NodePtr right{nullptr};
        std::weak_ptr<Node> parent;
    };


    std::shared_ptr<Node> root_;
};