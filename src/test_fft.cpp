#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "disfft.h"

using Node = SplittingTree::Node;
using NodePtr = SplittingTree::NodePtr;

bool CheckEqual(complex_t a, complex_t b) {
    return !NonZero(a - b);
}

TEST_CASE("Filters frequency simple 1") {
    auto root = std::make_shared<Node>();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();

    auto filter_a = Filter(a, 2);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(1), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(0), 0.));

    auto filter_b = Filter(b, 2);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(1), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(0), 1.));
}

TEST_CASE("Filters frequency simple 2") {
    {
        auto root = std::make_shared<Node>();
        auto a = root->MakeLeft();

        auto filter_a = Filter(a, 2);
        REQUIRE(CheckEqual(filter_a.FilterFrequency(1), 1));
    }
    {
        auto root = std::make_shared<Node>();
        auto b = root->MakeRight();

        auto filter_b = Filter(b, 2);
        REQUIRE(CheckEqual(filter_b.FilterFrequency(0), 1));
    }
}


TEST_CASE("Filters frequency simple 3") {
    auto root = std::make_shared<Node>();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();

    auto filter_a = Filter(a, 4);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(1), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(3), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(0), 0.));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(2), 0.));

    auto filter_b = Filter(b, 4);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(1), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(3), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(0), 1.));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(2), 1.));
}

TEST_CASE("Filters frequency simple 4") {
    auto root = std::make_shared<Node>();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();
    auto a1 = a->MakeLeft();
    auto b2 = b->MakeRight();

    auto filter_a = Filter(a1, 4);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(3), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(0), 0.));

    auto filter_b = Filter(b2, 4);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(3), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(0), 1.));
}

TEST_CASE("Filters frequency simple 5") {
    auto root = std::make_shared<Node>();
    auto a = root->MakeLeft();
    auto a1 = a->MakeLeft();
    auto a2 = a->MakeRight();

    auto filter_a = Filter(a1, 4);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(3), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(1), 0.));

    auto filter_b = Filter(a2, 4);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(3), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(1), 1.));
}

TEST_CASE("Filters time simple 1") {
    auto root = std::make_shared<Node>();
    auto a = root->MakeLeft();
    auto filter = Filter(a, 2);
    auto it = filter.FilterTime().find(0);
    REQUIRE(it != filter.FilterTime().end());
    REQUIRE(it->second == 1.);
    REQUIRE(filter.FilterTime().find(1) == filter.FilterTime().end());
}

TEST_CASE("Filters time simple 2") {
    auto root = std::make_shared<Node>();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();
    auto filter = Filter(a, 2);
    auto& filter_time = filter.FilterTime();
    auto it = filter_time.find(0);
    REQUIRE(it != filter_time.end());
    REQUIRE(it->second == 0.5);
    it = filter_time.find(1);
    REQUIRE(it != filter_time.end());
    REQUIRE(CheckEqual(it->second, -0.5));
}

TEST_CASE("Tree test remove") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();
    REQUIRE(root->left == a);
    REQUIRE(root->right == b);
    REQUIRE(a->parent == root.get());
    REQUIRE(b->parent == root.get());
    REQUIRE(root->parent == nullptr);
    REQUIRE(a->left == nullptr);
    REQUIRE(a->right == nullptr);
    REQUIRE(b->left == nullptr);
    REQUIRE(b->right == nullptr);

    tree.RemoveNode(a);
    REQUIRE(root->left == nullptr);
    REQUIRE(root->right == b);
    REQUIRE(b->parent == root.get());
    REQUIRE(root->parent == nullptr);
    REQUIRE(b->left == nullptr);
    REQUIRE(b->right == nullptr);

    tree.RemoveNode(b);
    REQUIRE(!tree.IsNonEmpty());
    REQUIRE(tree.GetRoot() == nullptr);
}

TEST_CASE("Tree get lightest 1") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    REQUIRE(root == tree.GetLightestNode());
    auto a = root->MakeLeft();
    REQUIRE(a == tree.GetLightestNode());
    auto b = root->MakeRight();
    auto node = tree.GetLightestNode();
    REQUIRE(b == node);
    tree.RemoveNode(b);
    node = tree.GetLightestNode();
    REQUIRE(node == a);
}

TEST_CASE("Tree get lightest 2") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a1 = root->MakeLeft();
    auto a2 = root->MakeRight();
    auto b1 = a1->MakeLeft();
    auto b2 = a1->MakeRight();
    auto c1 = b2->MakeLeft();
    auto c2 = b2->MakeRight();
    REQUIRE(a2 == tree.GetLightestNode());
    tree.RemoveNode(a2);
    REQUIRE(tree.GetLightestNode() == b1);
    tree.RemoveNode(b1);
    REQUIRE(c2 == tree.GetLightestNode());
    tree.RemoveNode(c2);
    REQUIRE(c1 == tree.GetLightestNode());
    tree.RemoveNode(c1);
    REQUIRE(!tree.IsNonEmpty());
}
