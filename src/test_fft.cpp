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