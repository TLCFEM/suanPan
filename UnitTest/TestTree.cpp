#include "CatchHeader.h"

#include <Toolbox/tree.h>

TEST_CASE("Quadtree constructs and accepts a single point", "[Utility.Tree]") {
    QuadTree<1> tree({{0.0, 0.0}, {10.0, 10.0}});
    REQUIRE_NOTHROW(tree.insert(Node2D{1.0, -1.0}));
}

TEST_CASE("Quadtree accepts points on boundaries", "[Utility.Tree]") {
    QuadTree<1> tree({{0.0, 0.0}, {2.0, 3.0}});

    REQUIRE_NOTHROW(tree.insert(Node2D{2.0, 0.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D{-2.0, 0.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D{0.0, 3.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D{0.0, -3.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D{2.0, 3.0}));
    REQUIRE_NOTHROW(tree.insert(Node2D{-2.0, -3.0}));
}

TEST_CASE("Quadtree handles many insertions across quadrants", "[Utility.Tree]") {
    QuadTree<2> tree({{0.0, 0.0}, {100.0, 100.0}});

    std::vector<Node2D> points{
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
    QuadTree<64> tree({{0.0, 0.0}, {1.0, 1.0}});

    for(unsigned i = 0; i < 32; ++i) REQUIRE_NOTHROW(tree.insert(Node2D{0.125, -0.25}));
}
