#include "catch.hpp"
#include <unordered_set>
#include <random>
#include <algorithm>
#include <set>

#include "utility_test.h"


TEST_CASE("Filters frequency simple 1") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();

    SignalInfo info(1, 2);

    auto filter_a = Filter(tree, a, info);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {1}}), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {0}}), 0.));

    auto filter_b = Filter(tree, b, info);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {1}}), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {0}}), 1.));
}

TEST_CASE("Filters frequency simple 2") {
    SignalInfo info(1, 2);
    {
        auto tree = SplittingTree();
        auto root = tree.GetRoot();
        auto a = root->MakeLeft();

        auto filter_a = Filter(tree, a, info);
        REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {1}}), 1));
    }
    {
        auto tree = SplittingTree();
        auto root = tree.GetRoot();
        auto b = root->MakeRight();

        auto filter_b = Filter(tree, b, info);
        REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {0}}), 1));
    }
}


TEST_CASE("Filters frequency simple 3") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();
    SignalInfo info(1, 4);

    auto filter_a = Filter(tree, a, info);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {1}}), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {3}}), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {0}}), 0.));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {2}}), 0.));

    auto filter_b = Filter(tree, b, info);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {1}}), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {3}}), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {0}}), 1.));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {2}}), 1.));
}

TEST_CASE("Filters frequency simple 4") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();
    auto a1 = a->MakeLeft();
    auto b2 = b->MakeRight();

    SignalInfo info(1, 4);

    auto filter_a = Filter(tree, a1, info);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {3}}), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {0}}), 0.));

    auto filter_b = Filter(tree, b2, info);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {3}}), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {0}}), 1.));
}

TEST_CASE("Filters frequency simple 5") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a = root->MakeLeft();
    auto a1 = a->MakeLeft();
    auto a2 = a->MakeRight();

    SignalInfo info(1, 4);

    auto filter_a = Filter(tree, a1, info);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {3}}), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(Key{info, {1}}), 0.));

    auto filter_b = Filter(tree, a2, info);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {3}}), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(Key{info, {1}}), 1.));
}

TEST_CASE("Filters frequency 2d 1") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a1 = root->MakeLeft();
    auto a2 = root->MakeRight();
    auto b1 = a2->MakeLeft();
    auto b2 = a2->MakeRight();
    auto c1 = b2->MakeLeft();
    auto c2 = b2->MakeRight();
    auto d1 = c1->MakeLeft();
    auto d2 = c1->MakeRight();

    SignalInfo info(2, 4);

    auto filter_d1 = Filter(tree, d1, info);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int value = (i == 0) && (j == 3);
            REQUIRE(CheckEqual(filter_d1.FilterFrequency(Key{info, {i, j}}), value));
        }
    }
}

TEST_CASE("Filters time simple 1") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a = root->MakeLeft();

    SignalInfo info(1, 2);

    auto filter = Filter(tree, a, info);
    auto it = filter.FilterValueAtTime(Key{info, {0}});
    REQUIRE(it == 1.);
    REQUIRE(filter.FilterValueAtTime(Key{info, {1}}) == 0.);
}

TEST_CASE("Filters time simple 2") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a = root->MakeLeft();
    auto b = root->MakeRight();
    SignalInfo info(1, 2);

    auto filter = Filter(tree, a, info);
    auto it = filter.FilterValueAtTime(Key{info, {0}});
    REQUIRE(it == 0.5);
    it = filter.FilterValueAtTime(Key{info, {1}});
    REQUIRE(CheckEqual(it, -0.5));
}

TEST_CASE("Filter full simple") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a1 = root->MakeLeft();
    auto a2 = root->MakeRight();
    auto b1 = a2->MakeLeft();
    auto b2 = a2->MakeRight();
    SignalInfo info(1, 4);

    auto filter = Filter(tree, b2, info);
    for (int i = 0; i < 4; ++i) {
        REQUIRE(CheckEqual(filter.FilterValueAtTime(Key{info, {i}}), 0.25));
    }
    REQUIRE(CheckEqual(filter.FilterFrequency(Key{info, {0}}), 1.));
    for (int i = 1; i < 4; ++i) {
        REQUIRE(CheckEqual(filter.FilterFrequency(Key{info, {i}}), 0.));
    }
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
    REQUIRE(tree.IsEmpty());
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
    REQUIRE(tree.IsEmpty());
}

std::pair<int, NodePtr> CustomGetLightest(NodePtr node) {
    if (!node->left && !node->right) {
        return {0, node};
    }
    std::vector<std::pair<int, NodePtr>> results;
    for (auto& child : {node->right, node->left}) {
        if (child) {
            results.push_back(CustomGetLightest(child));
        }
    }
    if (results.size() == 2) {
        if (results[0].first == results[1].first) {
            return {results[0].first + 1, results[0].second};
        } else {
            std::sort(results.begin(), results.end());
            return {results[0].first + 1, results[0].second};
        }
    } else {
        return *results.begin();
    }
}

TEST_CASE("Huge tree test 1") {
    auto tree = SplittingTree();
    std::unordered_set<NodePtr> leafs;
    leafs.insert(tree.GetRoot());
    const int iterations_num = 1e3;
    for (int i = 0; i < iterations_num; ++i) {
        auto r = rand() % leafs.size();
        auto it = leafs.begin();
        std::advance(it, r);
        auto random_leaf = *it;
        leafs.erase(random_leaf);
        auto res = random_leaf->AddChildren();
        leafs.insert(res.first);
        leafs.insert(res.second);
    }
    for (int i = 0; i < iterations_num + 1; ++i) {
        auto node = tree.GetLightestNode();
        int weight;
        NodePtr custom_node;
        std::tie(weight, custom_node) = CustomGetLightest(tree.GetRoot());
        REQUIRE(node.get() == custom_node.get());
        tree.RemoveNode(node);
    }
    REQUIRE(tree.IsEmpty());
}

TEST_CASE("DataSignal") {
    complex_t data[] = {0, 1, 2, 3, 4, 5, 6, 7};
    SignalInfo info(1, 8);

    DataSignal x(info, data);
    Key pos{info, 4};
    Key one{info, 1};
    int i = 4;
    for (int t = 0; t < 30; ++t) {
        REQUIRE(x.ValueAtTime(pos) == data[i]);
        i = (i - 1 + 8) % 8;
        pos.StoreDifference(pos, one);
    }
}

TEST_CASE("Key multidim") {
    SignalInfo info(3, 8);
    Key a(info, {4, 5, 3}), b(info, {7, 4, 7});
    Key buf(info);
    buf.StoreDifference(a, b);
    REQUIRE(buf == Key(info, {5, 1, 4}));
    REQUIRE(a * b == 4 * 7 + 5 * 4 + 7 * 3);
    REQUIRE(a.IncreaseAt(0, 3) == Key(info, {7, 5, 3}));
    REQUIRE(a.IncreaseAt(1, 3) == Key(info, {4, 0, 3}));
    REQUIRE(a.IncreaseAt(2, 3) == Key(info, {4, 5, 6}));
    REQUIRE(-a == Key(info, {4, 3, 5}));
}

TEST_CASE("ZeroTest 1") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a1 = root->MakeLeft();
    auto a2 = root->MakeRight();
    FrequencyMap chi{};
    SignalInfo info(1, 4);
    IndexGenerator delta(info, 321);
    TransformSettings settings;

    {
        // ifft([1, 0, 0, 0])
        complex_t data[] = {0.25 + 0.i, 0.25 - 0.i, 0.25 + 0.i, 0.25 + 0.i};
        DataSignal x(info, data);

        REQUIRE(!ZeroTest(x, chi, tree, a1, info, 1, delta, settings));
        REQUIRE(ZeroTest(x, chi, tree, a2, info, 1, delta, settings));
    }
    {
        // ifft([1, 0, 1, 0])
        complex_t data[] = {0.5+0.i, 0. -0.i, 0.5+0.i, 0. +0.i};
        DataSignal x(info, data);

        REQUIRE(!ZeroTest(x, chi, tree, a1, info, 2, delta, settings));
        REQUIRE(ZeroTest(x, chi, tree, a2, info, 2, delta, settings));
    }
    {
        // ifft([0, 1, 1, 0])
        complex_t data[] = {0.5 +0.i  , -0.25+0.25i,  0.  +0.i  , -0.25-0.25i};
        DataSignal x(info, data);

        REQUIRE(ZeroTest(x, chi, tree, a1, info, 2, delta, settings));
        REQUIRE(ZeroTest(x, chi, tree, a2, info, 2, delta, settings));
    }
}

TEST_CASE("ZeroTest 1.5") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a1 = root->MakeLeft();
    auto a2 = root->MakeRight();
    auto b1 = a2->MakeLeft();
    auto b2 = a2->MakeRight();
    FrequencyMap chi{};
    SignalInfo info(1, 4);
    IndexGenerator delta(info, 321);
    TransformSettings settings;

    {
        // ifft([0, 1, 1, 0])
        complex_t data[] = {0.5 +0.i  , -0.25+0.25i,  0.  +0.i  , -0.25-0.25i};
        DataSignal x(info, data);

        REQUIRE(ZeroTest(x, chi, tree, b1, info, 2, delta, settings));
        REQUIRE(!ZeroTest(x, chi, tree, b2, info, 2, delta, settings));
    }
}

TEST_CASE("ZeroTest 2") {
    auto tree = SplittingTree();
    auto root = tree.GetRoot();
    auto a = root->MakeRight();
    auto b1 = root->MakeLeft();
    auto b2 = b1->MakeLeft();
    auto b3 = b2->MakeRight();
    SignalInfo info(1, 8);
    IndexGenerator delta(info, 321);
    TransformSettings settings;

    {
        // ifft([0, 0, 0, 0, 0, 0, 0, 0])
        complex_t data[] = {0, 0, 0, 0, 0, 0, 0, 0};
        FrequencyMap chi{};
        DataSignal x(info, data);

        REQUIRE(!ZeroTest(x, chi, tree, a, info, 1, delta, settings));
        REQUIRE(!ZeroTest(x, chi, tree, b3, info, 1, delta, settings));
    }
    {
        // ifft([0, 0, 0, 1, 0, 0, 0, 0])
        complex_t data[] = {0.125     +0.i        , -0.08838835+0.08838835i,
                            0.        -0.125i     ,  0.08838835+0.08838835i,
                            -0.125     +0.i        ,  0.08838835-0.08838835i,
                            0.        +0.125i     , -0.08838835-0.08838835i };
        FrequencyMap chi{};
        DataSignal x(info, data);

        REQUIRE(!ZeroTest(x, chi, tree, a, info, 1, delta, settings));
        REQUIRE(ZeroTest(x, chi, tree, b3, info, 1, delta, settings));
    }
    {
        // ifft([0, 0, 1, 1, 1, 0, 0, 0])
        complex_t data[] = { 0.375     +0.i        , -0.21338835+0.21338835i,
                             0.        -0.125i     , -0.03661165-0.03661165i,
                             0.125     +0.i        , -0.03661165+0.03661165i,
                             0.        +0.125i     , -0.21338835-0.21338835i};
        FrequencyMap chi{};
        DataSignal x(info, data);

        REQUIRE(ZeroTest(x, chi, tree, a, info, 3, delta, settings));
        REQUIRE(ZeroTest(x, chi, tree, b3, info, 3, delta, settings));
    }
    {
        // ifft([0, 0, 0, 1, 0, 0, 0, 1])
        complex_t data[] = { 0.25+0.i  ,  0.  -0.i  ,  0.  -0.25i,  0.  -0.i  , -0.25+0.i  ,
                             0.  +0.i  ,  0.  +0.25i,  0.  +0.i  };
        FrequencyMap chi;
        chi[Key{info, {7}}] = 1.;
        DataSignal x(info, data);

        REQUIRE(!ZeroTest(x, chi, tree, a, info, 2, delta, settings));
        REQUIRE(ZeroTest(x, chi, tree, b3, info, 2, delta, settings));
    }
    {
        // ifft([0, 0, 0, 2, 0, -1, 0, 3])
        complex_t data[] = { 0.5      +0.00000000e+00i,  0.1767767-1.38777878e-17i,
                             0.       -7.50000000e-01i, -0.1767767-1.38777878e-17i,
                             -0.5      +0.00000000e+00i, -0.1767767+1.38777878e-17i,
                             0.       +7.50000000e-01i,  0.1767767+1.38777878e-17i};
        FrequencyMap chi;
        chi[Key{info, {7}}] = 3.;
        chi[Key{info, {5}}] = -1.;
        DataSignal x(info, data);

        REQUIRE(!ZeroTest(x, chi, tree, a, info, 3, delta, settings));
        REQUIRE(ZeroTest(x, chi, tree, b3, info, 3, delta, settings));
    }
    {
        // ifft([0, 0, 0, 0, 0, -1, 0, 3])
        complex_t data[] = { 0.25               +0.i                 ,
                             0.35355339059327373 - 0.17677669529663692i,
                             0.                 -0.5i                ,
                             -0.35355339059327373 - 0.17677669529663692i,
                             -0.25               +0.i                 ,
                             -0.35355339059327373 + 0.17677669529663692i,
                             0.                 +0.5i                ,
                             0.35355339059327373 + 0.17677669529663692i};
        FrequencyMap chi;
        chi[Key{info, {7}}] = 3.;
        chi[Key{info, {5}}] = -1.;
        DataSignal x(info, data);

        REQUIRE(!ZeroTest(x, chi, tree, a, info, 2, delta, settings));
        REQUIRE(!ZeroTest(x, chi, tree, b3, info, 2, delta, settings));
    }
}

TEST_CASE("SparseFFT 1") {
    SignalInfo info{1, 4};
    REQUIRE(RunSFFT(info, 0,
            {0, 0, 0, 0},
            {0, 0, 0, 0}
    ));
    REQUIRE(RunSFFT(info, 1,
            { 0.25+0.i  ,  0.  +0.25i, -0.25+0.i  ,  0.  -0.25i},
            {0, 1, 0, 0}
    ));
    REQUIRE(RunSFFT(info, 2,
            {0.5 +0.i  , -0.25+0.25i,  0.  +0.i  , -0.25-0.25i},
            {0, 1, 1, 0}
    ));
    REQUIRE(RunSFFT(info, 2,
            {0.5 +0.i  , -0.25-0.25i,  0.  +0.i  , -0.25+0.25i},
            {0, 0, 1, 1}
    ));
    REQUIRE(RunSFFT(info, 3,
            {0.75+0.i  , 0.  -0.25i, 0.25+0.i  , 0.  +0.25i},
            {1, 0, 1, 1}
    ));
    REQUIRE(RunSFFT(info, 4,
            {1, 0, 0, 0},
            {1, 1, 1, 1}
    ));
}

TEST_CASE("SparseFFT 2") {
    REQUIRE(RunSFFT({1, 8}, 3,
                    { 1.25               +0.i                 ,
                      -0.8169417382415922 -0.44194173824159216i,
                      0.625              +0.625i              ,
                      0.06694173824159222-0.44194173824159216i,
                      0.                 +0.i                 ,
                      0.06694173824159222+0.44194173824159216i,
                      0.625              -0.625i              ,
                      -0.8169417382415922 +0.44194173824159216i},
                    {1, 0, 0, 0, 4, 5, 0, 0}
    ));
}
