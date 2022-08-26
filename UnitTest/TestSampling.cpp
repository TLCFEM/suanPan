#include <Toolbox/resampling.h>
#include "CatchHeader.h"

TEST_CASE("GCD", "[Utility.Sampling]") {
    REQUIRE(gcd(20, 10) == 10);
    REQUIRE(gcd(14, 24) == 2);
    REQUIRE(gcd(56, 98) == 14);
    REQUIRE(gcd(51248431, 3411571202) == 13);
}
