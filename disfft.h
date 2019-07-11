#pragma once

#define _USE_MATH_DEFINES
#define PI M_PI
#define I complex_t(0, 1)

#include "complex"

#include "assert.h"
#include "memory"
#include <vector>
#include <algorithm>

using complex_t  = std::complex<double>;

class SplittingTree {
public:
    struct Node;

    SplittingTree(): root_(std::make_shared<Node>()) {
    }

    using NodePtr = std::shared_ptr<Node>;

    struct Node {
        int length{0};
        uint64_t label{0};

        NodePtr left{nullptr};
        NodePtr right{nullptr};
//        std::weak_ptr<Node> parent;
        Node* parent{nullptr};

        Node(int length, uint64_t label) : length(length), label(label) {}

        Node CalcParent() const {
            assert(length > 0);
            return {length - 1, label & ((1ull << static_cast<uint64_t>(length - 1)) - 1)};
        }

        Node CalcLeft() const {
            assert(length > 0);
            return {length + 1, label + (1ull << static_cast<uint64_t>(length + 1))};
        }

        Node CalcRight() const {
            assert(length > 0);
            return {length + 1, label};
        }

        bool HasBothChild() const {
            return left && right;
        }


        std::vector<bool> GetRootPathMask () const {
            std::vector<bool> mask;
//            for (auto node = this; node != nullptr; node = node->parent.lock().get()) {
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