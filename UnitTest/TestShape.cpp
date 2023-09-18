#include <Toolbox/shape.h>
#include "CatchHeader.h"

TEST_CASE("Compute Area By Shoelace", "[Utility.Shape]") {
    const mat C{{3, 4}, {5, 6}, {9, 5}, {12, 8}, {5, 11}};

    REQUIRE(Approx(30) == area::shoelace(C));
}

TEST_CASE("Small Sparse", "[Utility.Shape]") {
    for(auto I = 2; I <= 20; I += 2) {
        const auto B = sprandu(6, 6 * I, .3);
        const auto D = mat(6, 6, fill::randn);
        const auto BB = mat(B);

        BENCHMARK("Sparse BTDB") {
            mat A = B.t() * D * B;
            return A;
        };

        BENCHMARK("Dense BTDB") {
            mat A = BB.t() * D * BB;
            return A;
        };
    }
}
