#include "CatchHeader.h"

#include <Toolbox/container.h>

TEST_CASE("Set to uvec", "[Utility.Container]") {
    const suanpan::set<int> s{3, 1, 2};
    const auto v = to_uvec(s);
    REQUIRE(3 == v.n_elem);
    REQUIRE(1 == v(0));
    REQUIRE(2 == v(1));
    REQUIRE(3 == v(2));
}

TEST_CASE("Unordered Set to uvec", "[Utility.Container]") {
    const suanpan::unordered_set<int> s{5, 3, 8};
    const auto v = to_uvec(s);
    REQUIRE(3 == v.n_elem);
    // order of unordered set is unspecified; only membership matters
    uvec vals = sort(v);
    REQUIRE(3 == vals(0));
    REQUIRE(5 == vals(1));
    REQUIRE(8 == vals(2));
}

TEST_CASE("Container Aliases", "[Utility.Container]") {
    suanpan::vector<int> v{1, 2, 3};
    REQUIRE(3 == v.size());
    v.emplace_back(4);
    REQUIRE(4 == v.size());
    REQUIRE(4 == v.at(3));

    suanpan::unordered_map<int, int> um;
    um[10] = 100;
    REQUIRE(100 == um[10]);

    suanpan::set<int> s{1, 2};
    s.insert(3);
    REQUIRE(3 == s.size());
    REQUIRE(1 == s.count(2));

    suanpan::unordered_set<int> us{1, 2, 3};
    REQUIRE(3 == us.size());

    suanpan::graph<int> g(3);
    g[0].insert(1);
    g[1].insert(0);
    REQUIRE(1 == g[0].size());
}
