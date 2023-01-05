#include <Toolbox/shape.h>
#include "CatchHeader.h"

TEST_CASE("Compute Area By Shoelace", "[Utility.Shape]") {
    const mat C{{3, 4}, {5, 6}, {9, 5}, {12, 8}, {5, 11}};

    REQUIRE(Approx(30) == area::shoelace(C));
}
