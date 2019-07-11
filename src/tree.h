#pragma once

#include "assert.h"
#include "memory"
#include "vector"
#include "arithmetics.h"

class SplittingTree {
public:
    struct Node;

    SplittingTree(): root_(std::make_shared<Node>()) {
    }

    using NodePtr = std::shared_ptr<Node>;

    struct Node {
        int level{0};
        uint64_t label{0};

        NodePtr left{nullptr};
        NodePtr right{nullptr};
        Node* parent{nullptr};

        Node(int length, uint64_t label) : level(length), label(label) {}
        Node() {}

        Node CalcParent() const {
            assert(level > 0);
            return {level - 1, label & ((1ull << static_cast<uint64_t>(level - 1)) - 1)};
        }

        Node CalcLeft() const {
            assert(level > 0);
            return {level + 1, label + (1ull << static_cast<uint64_t>(level + 1))};
        }

        Node CalcRight() const {
            assert(level > 0);
            return {level + 1, label};
        }

        bool HasBothChild() const {
            return left && right;
        }


        std::vector<bool> GetRootPathMask() const {
            std::vector<bool> mask;
            for (auto node = this; node != nullptr; node = node->parent) {
                mask.push_back(node->HasBothChild());
            }
            return mask;
        }
    };

    NodePtr GetRoot() const {
        return root_;
    }

private:
    NodePtr root_;
};

