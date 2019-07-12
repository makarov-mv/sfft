#pragma once

#include "assert.h"
#include "memory"
#include "vector"
#include "algorithm"
#include "arithmetics.h"

class SplittingTree {
public:
    struct Node;

    SplittingTree(): root_(std::make_shared<Node>()) {
    }

    using NodePtr = std::shared_ptr<Node>;

    struct Node {
        int level{0};
        int64_t label{0};

        NodePtr left{nullptr};
        NodePtr right{nullptr};
        Node* parent{nullptr};

        Node(int length, uint64_t label) : level(length), label(label) {}
        Node() {}

        NodePtr MakeLeft() const {
            assert(level > 0);
            return std::make_shared<Node>(level + 1, label + (1ull << static_cast<uint64_t>(level + 1)));
        }

        NodePtr MakeRight() const {
            assert(level > 0);
            return std::make_shared<Node>(level + 1, label);
        }

        bool HasBothChild() const {
            return left && right;
        }

        void AddChildren() {
            left = MakeLeft();
            right = MakeRight();
        }


        std::vector<bool> GetRootPathMask() const {
            std::vector<bool> mask;
            for (auto node = this; node != nullptr; node = node->parent) {
                mask.push_back(node->HasBothChild());
            }
            std::reverse(mask.begin(), mask.end());
            return mask;
        }
    };

    NodePtr GetRoot() const {
        return root_;
    }

    NodePtr GetLightestNode() const {
        return getLightest(root_);
    }

    void RemoveNode(const NodePtr& node) {
        Node* current = node.get();
        while (current && !current->HasBothChild()) {
            current = current->parent;
        }
        if (!current) {
            root_ = nullptr;
        }
        if (current->left->parent == current) {
            current->left = nullptr;
        } else {
            current->right = nullptr;
        }
    }

    bool IsNonEmpty() const {
        return root_ != nullptr;
    }

private:

    NodePtr getLightest(const NodePtr& node) const {
        if (!node->left && !node->right) {
            return node;
        }
        NodePtr lightest(nullptr);
        if (node->left) {
            lightest = getLightest(node->left);
        }
        if (node->right) {
            if (lightest) {
                auto right_lightest = getLightest(node->right);
                lightest = lightest->level < right_lightest->level ? lightest : right_lightest;
            } else {
                lightest = getLightest(node->right);
            }
        }
        return lightest;
    }

    NodePtr root_;
};

