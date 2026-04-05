#include "CatchHeader.h"

#include <Toolbox/tree/pquadtree.hpp>
#include <Toolbox/tree/quadtree.hpp>
#include <random>

TEST_CASE("Quadtree constructs and accepts a single point", "[Utility.Tree]") {
    QuadTree tree({{0.0, 0.0}, {10.0, 10.0}});
    REQUIRE_NOTHROW(tree.insert(Node2D<>{1.0, -1.0}));
}

TEST_CASE("Quadtree accepts points on boundaries", "[Utility.Tree]") {
    QuadTree tree({{0.0, 0.0}, {2.0, 3.0}});

    REQUIRE_NOTHROW(tree.insert(Node2D<>{2.0, 0.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D<>{-2.0, 0.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D<>{0.0, 3.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D<>{0.0, -3.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D<>{2.0, 3.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D<>{-2.0, -3.0}));
}

TEST_CASE("Quadtree handles many insertions across quadrants", "[Utility.Tree]") {
    QuadTree<double, 2> tree({{0.0, 0.0}, {100.0, 100.0}});

    std::vector<Node2D<>> points{
        {-50.0, -50.0},
        {-25.0, -25.0},
        {-10.0, -60.0},
        {-70.0, -15.0},
        {50.0, -50.0},
        {25.0, -25.0},
        {10.0, -60.0},
        {70.0, -15.0},
        {50.0, 50.0},
        {25.0, 25.0},
        {10.0, 60.0},
        {70.0, 15.0},
        {-50.0, 50.0},
        {-25.0, 25.0},
        {-10.0, 60.0},
        {-70.0, 15.0},
    };

    for(auto&& p : points) REQUIRE_NOTHROW(tree.insert(std::move(p)));
}

TEST_CASE("Quadtree handles repeated insertions without splitting when bucket is large", "[Utility.Tree]") {
    QuadTree<double, 64> tree({{0.0, 0.0}, {1.0, 1.0}});

    for(unsigned i = 0; i < 32; ++i) REQUIRE_NOTHROW(tree.insert(Node2D<>{0.125, -0.25}));
}

TEST_CASE("Quadtree benchmark: insert random nodes", "[Utility.Tree]") {
    static constexpr auto size = 1000.;
    static constexpr auto n = 30'000;
    std::mt19937 gen(42);
    std::uniform_real_distribution dis(-size, size);

    BENCHMARK("Insert random nodes") {
        QuadTree<double, 2> tree({{0.0, 0.0}, {size, size}});
        for(auto i = 0; i < n; ++i) tree.insert(Node2D<>{dis(gen), dis(gen)});
        return tree;
    };
}

TEST_CASE("PQuadtree benchmark: insert random nodes", "[Utility.Tree]") {
    static constexpr auto size = 1000.;
    static constexpr auto n = 30'000;
    std::mt19937 gen(42);
    std::uniform_real_distribution dis(-size, size);

    using tree_t = PQuadTree<double, 2>;
    using node_t = tree_t::node_t;

    std::vector<node_t> points;
    points.reserve(n);
    for(auto i = 0; i < n; ++i) points.push_back(node_t{dis(gen), dis(gen)});

    std::vector<tree_t::node_ptr_t> ptrs;
    ptrs.reserve(points.size());
    for(auto&& p : points) ptrs.push_back(&p);

    BENCHMARK("Insert random nodes") {
        tree_t tree({{0.0, 0.0}, {size, size}});
        tree.insert(std::move(ptrs));
        tree.overlap({{dis(gen), dis(gen)}, {1, 1}});
        return tree;
    };
}
