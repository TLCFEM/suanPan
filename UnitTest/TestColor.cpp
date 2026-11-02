#include "CatchHeader.h"

#include <Toolbox/sort_color.hpp>

namespace {
    /**
     * @brief Build a graph of N nodes where every node gets a self-loop and the given
     *        inter-node edges are added symmetrically.
     *
     * Note: the coloring algorithms require every node to possess at least one
     * neighbour (the self-loop) in order to terminate, mirroring the convention used
     * by the rest of the test suite.
     */
    auto make_graph(const uword N) {
        suanpan::graph<unsigned> g(N);
        for(uword I = 0; I < N; ++I) g[I].insert(static_cast<unsigned>(I)); // self loop
        return g;
    }

    void connect(suanpan::graph<unsigned>& g, const uword A, const uword B) {
        g[A].insert(static_cast<unsigned>(B));
        g[B].insert(static_cast<unsigned>(A));
    }

    /**
     * @brief Validate that a coloring is valid for the given graph.
     *
     * A valid proper graph coloring must:
     *  - assign every node exactly once;
     *  - never assign two adjacent (connected, distinct) nodes the same color.
     */
    void validate_coloring(const suanpan::graph<unsigned>& node_register, const std::vector<std::vector<unsigned>>& color_map) {
        std::vector seen(node_register.size(), false);

        for(const auto& color : color_map) {
            REQUIRE_FALSE(color.empty());
            for(const auto node : color) {
                REQUIRE(node < node_register.size());
                REQUIRE_FALSE(seen[node]);
                seen[node] = true;
            }
        }

        for(const auto& color : color_map)
            for(auto I = 0u; I < color.size(); ++I)
                for(auto J = I + 1u; J < color.size(); ++J)
                    // nodes in a color are all distinct (see `seen` above), so a self-loop is impossible
                    REQUIRE_FALSE(node_register[color[I]].contains(color[J]));

        for(auto flag : seen) REQUIRE(flag);
    }
} // namespace

TEST_CASE("Coloring Empty Graph", "[Utility.Color]") {
    const suanpan::graph<unsigned> empty{};
    REQUIRE(sort_color_wp(empty).empty());
    REQUIRE(sort_color_mis(empty).empty());
}

TEST_CASE("Coloring Single Node Graph", "[Utility.Color]") {
    const auto g = make_graph(1);
    const auto wp = sort_color_wp(g);
    REQUIRE(1 == wp.size());
    REQUIRE(1 == wp[0].size());
    REQUIRE(0u == wp[0][0]);

    const auto mis = sort_color_mis(g);
    REQUIRE(1 == mis.size());
    REQUIRE(1 == mis[0].size());
}

TEST_CASE("Coloring Complete Graph", "[Utility.Color]") {
    // a complete graph (every pair of distinct nodes connected) needs N colors
    constexpr uword N = 5;
    auto g = make_graph(N);
    for(uword I = 0; I < N; ++I)
        for(uword J = I + 1; J < N; ++J) connect(g, I, J);

    const auto wp = sort_color_wp(g);
    REQUIRE(N == wp.size());
    validate_coloring(g, wp);

    const auto mis = sort_color_mis(g);
    REQUIRE(N == mis.size());
    validate_coloring(g, mis);
}

TEST_CASE("Coloring Nodes With Only Self Loops", "[Utility.Color]") {
    // nodes with only self-loops are effectively independent: one color suffices
    constexpr uword N = 7;
    const auto isolated = make_graph(N);

    const auto wp = sort_color_wp(isolated);
    REQUIRE(1 == wp.size());
    REQUIRE(N == wp[0].size());
    validate_coloring(isolated, wp);

    const auto mis = sort_color_mis(isolated);
    REQUIRE(1 == mis.size());
    REQUIRE(N == mis[0].size());
    validate_coloring(isolated, mis);
}

TEST_CASE("Coloring Bipartite Graph", "[Utility.Color]") {
    // a complete bipartite graph - greedy coloring must still be valid
    constexpr uword N = 6;
    auto g = make_graph(N);
    for(uword I = 0; I < N / 2; ++I)
        for(uword J = N / 2; J < N; ++J) connect(g, I, J);

    const auto wp = sort_color_wp(g);
    validate_coloring(g, wp);

    const auto mis = sort_color_mis(g);
    validate_coloring(g, mis);
}

TEST_CASE("Coloring Random Graphs", "[Utility.Color]") {
    for(auto trial = 0; trial < 20; ++trial) {
        constexpr uword N = 40;
        auto g = make_graph(N);
        for(uword I = 0; I < N; ++I)
            for(uword J = I + 1; J < N; ++J)
                if(randu() < .1) connect(g, I, J);

        const auto wp = sort_color_wp(g);
        validate_coloring(g, wp);
        const auto mis = sort_color_mis(g);
        validate_coloring(g, mis);
    }
}
