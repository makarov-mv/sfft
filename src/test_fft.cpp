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

    auto filter_a = Filter(a, 3);
    REQUIRE(CheckEqual(filter_a.FilterFrequency(1), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(3), 1));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(0), 0.));
    REQUIRE(CheckEqual(filter_a.FilterFrequency(2), 0.));

    auto filter_b = Filter(b, 2);
    REQUIRE(CheckEqual(filter_b.FilterFrequency(1), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(3), 0));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(0), 1.));
    REQUIRE(CheckEqual(filter_b.FilterFrequency(2), 1.));
}
