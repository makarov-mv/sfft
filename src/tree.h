#pragma once

#include "assert.h"
#include "memory"
#include "vector"
#include "algorithm"
#include "arithmetics.h"
#include "tuple"

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

        Node(int length, int64_t label) : level(length), label(label) {}
        Node() {}

        NodePtr MakeLeft() {
            assert(level >= 0);
            auto son = std::make_shared<Node>(level + 1, label + (1ll << static_cast<int64_t>(level)));
            son->parent = this;
            this->left = son;
            return son;
        }

        NodePtr MakeRight() {
            assert(level >= 0);
            auto son = std::make_shared<Node>(level + 1, label);
            son->parent = this;
            this->right = son;
            return son;
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
        return getLightest(root_).first;
    }

    void RemoveNode(const NodePtr& node) {
        Node* current = node.get();
        Node* son = nullptr;
        while (current && !current->HasBothChild()) {
            son = current;
            current = current->parent;
        }
        if (!current) {
            root_.reset();
            return;
        }
        if (current->left.get() == son) {
            current->left.reset();
        } else {
            current->right.reset();
        }
    }

    bool IsNonEmpty() const {
        return root_ != nullptr;
    }

private:
    std::pair<NodePtr, int> getLightest(const NodePtr& node) const {
        if (!node->left && !node->right) {
            return {node, 0};
        }
        NodePtr lightest(nullptr);
        int weight;
        if (node->left) {
            std::tie(lightest, weight) = getLightest(node->left);
        }
        if (node->right) {
            if (lightest) {
                auto [right_lightest, right_weight] = getLightest(node->right);
                lightest = weight < right_weight ? lightest : right_lightest;
                return {lightest, std::min(weight, right_weight) + 1};
            } else {
                std::tie(lightest, weight) = getLightest(node->right);
                return {lightest, weight};
            }
        }
        return {lightest, weight};
    }

    NodePtr root_;
};

