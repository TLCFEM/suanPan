#include "CatchHeader.h"

#include <Toolbox/utility.h>

TEST_CASE("String Compare", "[Utility.String]") {
    const std::string A = "logic";
    const std::string B = "LoGicAnD";
    const std::string C = "LoGic";

    REQUIRE(!is_equal(A, B));
    REQUIRE(is_equal(A, B.substr(0, 5)));
    REQUIRE(is_equal(A, C));
}

TEST_CASE("String View Compare", "[Utility.String]") {
    constexpr std::string_view SOURCE = "LoGicAndMore";

    REQUIRE(is_equal(SOURCE.substr(0, 5), "logic"));
    REQUIRE(!is_equal(SOURCE.substr(0, 6), "logic"));
    REQUIRE(!is_equal(SOURCE.substr(1, 5), "logic"));

    const std::string WITH_NULL{"Lo\0gic", 6};

    REQUIRE(is_equal(std::string_view{WITH_NULL}, std::string_view{"Lo\0gic", 6}));
    REQUIRE(!is_equal(std::string_view{WITH_NULL}, std::string_view{"Lo\0gid", 6}));
    REQUIRE(!is_equal(std::string_view{WITH_NULL}, std::string_view{"Lo\0gicX", 7}));
}
