#include <Toolbox/utility.h>
#include "CatchHeader.h"

TEST_CASE("String Compare", "[Utility.String]") {
    const std::string A = "logic";
    const std::string B = "LoGicAnD";
    const std::string C = "LoGic";

    REQUIRE(!is_equal(A, B));
    REQUIRE(is_equal(A, B.substr(0, 5)));
    REQUIRE(is_equal(A, C));
}
