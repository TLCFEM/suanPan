#include "CatchHeader.h"

#include <Toolbox/utility.h>

TEST_CASE("Binomial Compute Basic Function", "[Utility.Binomial]") {
    REQUIRE(165 == suanpan::binomial(11, 3));
    REQUIRE(6435 == suanpan::binomial(15, 7));
    REQUIRE(455 == suanpan::binomial(15, 3));
    REQUIRE(126 == suanpan::binomial(9, 4));
    REQUIRE(220 == suanpan::binomial(12, 3));
}

TEST_CASE("Sign", "[Utility.Sign]") {
    REQUIRE(Approx(1) == suanpan::sign(std::numeric_limits<double>::epsilon()));
    REQUIRE(Approx(-1) == suanpan::sign(-std::numeric_limits<double>::epsilon()));
    REQUIRE(Approx(0) == suanpan::sign(0));
}

TEST_CASE("Matrix Allocation", "[Utility.Matrix]") {
    BENCHMARK("Static Size 20") {
        mat::fixed<20, 20> A(fill::randn);
        A(10, 10) = 1.;
        REQUIRE(A.n_elem == 400);
        return A;
    };

    BENCHMARK("Dynamic Size 20") {
        mat A(20, 20, fill::randn);
        A(10, 10) = 1.;
        REQUIRE(A.n_elem == 400);
        return A;
    };
}

TEST_CASE("Color Print", "[Utility.Print]") {
    suanpan_info("TEST.\n");
    suanpan_highlight("TEST.\n");
    suanpan_debug("TEST.\n");
    suanpan_warning("TEST.\n");
    suanpan_error("TEST.\n");
    suanpan_fatal("TEST.\n");
    suanpan_info("TEST.\n", vec{1, 2, 3});
}
