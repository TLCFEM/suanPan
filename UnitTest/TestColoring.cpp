#include <Toolbox/sort_color.h>
#include "CatchHeader.h"

auto tiny_graph() {
    suanpan_register graph;
    graph.emplace_back(std::initializer_list<unsigned>{4, 2, 1, 0});
    graph.emplace_back(std::initializer_list<unsigned>{0, 2, 3, 1});
    graph.emplace_back(std::initializer_list<unsigned>{4, 3, 1, 0, 2});
    graph.emplace_back(std::initializer_list<unsigned>{1, 2, 5, 6, 3});
    graph.emplace_back(std::initializer_list<unsigned>{0, 2, 5, 4});
    graph.emplace_back(std::initializer_list<unsigned>{4, 3, 6, 5});
    graph.emplace_back(std::initializer_list<unsigned>{5, 3, 6});
    return graph;
}

auto small_graph(const unsigned long long N, const double D) {
    suanpan_register graph(N);

    for(auto I = 0u; I < N; ++I) graph[I].insert(I);

    sp_mat A = sprandu(N, N, D);

    for(auto I = A.begin(); I != A.end(); ++I) {
        graph[I.col()].insert(static_cast<unsigned>(I.row()));
        graph[I.row()].insert(static_cast<unsigned>(I.col()));
    }

    return graph;
}

auto color_graph(const unsigned long long N, const double D, vector<vector<unsigned>> (&algorithm)(const suanpan_register&)) {
    auto graph = small_graph(N, D);

    const auto color_map = algorithm(graph);

    auto number_element = 0llu;
    for(auto& color : color_map) {
        number_element += color.size();
        for(auto I = 0llu; I < color.size(); ++I) {
            const auto& target_list = graph[color[I]];
            for(auto J = I + 1llu; J < color.size(); ++J)
                REQUIRE((target_list.find(color[J]) == target_list.end()));
        }
    }

    REQUIRE(number_element == graph.size());
}

TEST_CASE("Coloring Basic", "[Utility.Coloring]") {
    auto graph = tiny_graph();

    auto color_map = sort_color_wp(graph);

    auto number_element = 0llu;
    for(auto& color : color_map) {
        number_element += color.size();
        for(auto I = 0llu; I < color.size(); ++I) {
            const auto& target_list = graph[color[I]];
            for(auto J = I + 1llu; J < color.size(); ++J)
                REQUIRE((target_list.find(color[J]) == target_list.end()));
        }
    }

    REQUIRE(number_element == graph.size());

    color_map = sort_color_mis(graph);

    number_element = 0llu;
    for(auto& color : color_map) {
        number_element += color.size();
        for(auto I = 0llu; I < color.size(); ++I) {
            const auto& target_list = graph[color[I]];
            for(auto J = I + 1llu; J < color.size(); ++J)
                REQUIRE((target_list.find(color[J]) == target_list.end()));
        }
    }

    REQUIRE(number_element == graph.size());
}

TEST_CASE("Coloring WP", "[Utility.Coloring]") {
    for(auto N = 16; N < 512; N *= 2)
        for(auto D = 1; D < 8; D *= 2)
            BENCHMARK(string("Graph Size " + std::to_string(N) + " Density " + std::to_string(D)).c_str()) { color_graph(N, D * 1E-2, sort_color_wp); };
}

TEST_CASE("Coloring MIS", "[Utility.Coloring]") {
    for(auto N = 16; N < 512; N *= 2)
        for(auto D = 1; D < 8; D *= 2)
            BENCHMARK(string("Graph Size " + std::to_string(N) + " Density " + std::to_string(D)).c_str()) { color_graph(N, D * 1E-2, sort_color_mis); };
}
