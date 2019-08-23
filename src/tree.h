#pragma once

#include "assert.h"
#include "memory"
#include "vector"
#include "algorithm"
#include "arithmetics.h"
#include "tuple"
#include "iostream"

class SplittingTree {
public:
    struct Node;
    using NodePtr = std::shared_ptr<Node>;

    explicit SplittingTree(const SplittingTree* parent_tree, const NodePtr& parent_node): parent_tree_(parent_tree),
        parent_node_(parent_node), root_(std::make_shared<Node>()) {
        if (parent_tree) {
            root_->label = parent_node_->label;
            root_->level = parent_node_->level;
        }
    }

    SplittingTree(): parent_tree_(nullptr), parent_node_(nullptr), root_(std::make_shared<Node>()) {}

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

        std::pair<NodePtr, NodePtr> AddChildren() {
            return {MakeLeft(), MakeRight()};
        }

        bool Removed() const {
            return this == parent;
        }
    };

    std::vector<int> GetRootPath(const NodePtr& v) const {
        std::vector<int> path;
        ReversedRootPath(v, path);
        return path;
    }

    NodePtr GetRoot() const {
        return root_;
    }

    NodePtr GetLightestNode() const {
        return getLightest(root_).first;
    }

    void RemoveNode(NodePtr node) {
        Node* current = node.get();
        Node* son = nullptr;
        while (current && !current->HasBothChild()) {
            son = current;
            current = current->parent;
        }
        if (!current) {
            node->parent = node.get();
            DeleteNode(root_);
            return;
        }
        NodePtr other_son;
        if (current->left.get() == son) {
            DeleteNode(current->left);
            other_son = current->right;
        } else {
            DeleteNode(current->right);
            other_son = current->left;
        }
        if (current != root_.get()) {
            auto parent = current->parent;
            other_son->parent = parent;
            if (parent->left.get() == current) {
                parent->left = other_son;
            } else {
                parent->right = other_son;
            }
        }
        node->parent = node.get();
    }

    bool IsEmpty() const {
        return root_ == nullptr;
    }

    int LeavesCount() const {
        return root_ ? SubtreeLeavesCount(root_) : 0;
    }

    int SubtreeLeavesCount(const NodePtr& node) const {
        if (!node->left && !node->right) {
            return 1;
        }
        int res = 0;
        if (node->left) {
            res += SubtreeLeavesCount(node->left);
        }
        if (node->right) {
            res += SubtreeLeavesCount(node->right);
        }
        return res;
    }

private:
    void ReversedRootPath(const NodePtr& v, std::vector<int>& path) const {
        for (auto node = v.get(); node != nullptr; node = node->parent) {
            if (node->HasBothChild()) {
                path.push_back(node->level);
            }
        }
        if (parent_tree_) {
            parent_tree_->ReversedRootPath(parent_node_, path);
        }
    }

    void DeleteNode(NodePtr& node) const {
        node->parent = node.get();
        node.reset();
    }

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
            }
        }
        return {lightest, weight + node->HasBothChild()};
    }
    const SplittingTree* parent_tree_;
    NodePtr parent_node_;
    NodePtr root_;
};

void PrintNodeAsDot(const SplittingTree::NodePtr& node) {
    if (!node) {
        return;
    }
    std::cout << "n" << node.get() << "[label=\"level: " << node->level << "\nlabel: " << node->label << "\"];\n";
    if (node->left) {
        std::cout << "n" << node.get() << " -> n" << node->left.get() <<";\n";
    }
    if (node->right) {
        std::cout << "n" << node.get() << " -> n" << node->right.get() <<";\n";
    }
    PrintNodeAsDot(node->left);
    PrintNodeAsDot(node->right);
}

void PrintTreeAsDot(const SplittingTree& tree, const std::string& name) {
    std::cout << "digraph " << name << " {\n";
    PrintNodeAsDot(tree.GetRoot());
    std::cout << "}\n\n";
}