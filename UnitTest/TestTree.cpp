#include "CatchHeader.h"

#include <Toolbox/tree/quadtree.hpp>
#include <random>

TEST_CASE("Quadtree Benchmark", "[Utility.Tree]") {
    static constexpr auto size = 1000.;
    static constexpr auto n = 30'000;
    std::mt19937 gen(42);
    std::uniform_real_distribution dis(-size, size);

    using tree_t = QuadTree<double, 2>;
    using node_t = tree_t::node_t;

    std::vector<node_t> points;
    points.reserve(n);
    for(auto i = 0; i < n; ++i) points.push_back(node_t{dis(gen), dis(gen)});

    BENCHMARK("Random Insertion") {
        std::vector<tree_t::node_ptr_t> ptrs;
        ptrs.reserve(points.size());
        for(auto&& p : points) ptrs.push_back(&p);
        tree_t tree({{0.0, 0.0}, {size, size}});
        tree.insert(std::move(ptrs));
        tree.overlap({{dis(gen), dis(gen)}, {1, 1}});
        return tree;
    };
}
