#include <Toolbox/tensorToolbox.h>
#include "CatchHeader.h"

TEST_CASE("Invariant", "[Utility.Tensor]") {
    const vec A{50, -20, 10, 30, -10, 20};

    REQUIRE(Approx(40) == tensor::stress::invariant1(A));
    REQUIRE(Approx(2100) == tensor::stress::invariant2(A));
    REQUIRE(Approx(-28000) == tensor::stress::invariant3(A));

    const auto B = eig_sym(tensor::stress::to_tensor(A)).eval();

    REQUIRE(Approx(40) == tensor::stress::invariant1(B));
    REQUIRE(Approx(2100) == tensor::stress::invariant2(B));
    REQUIRE(Approx(-28000) == tensor::stress::invariant3(B));
}

TEST_CASE("Lode Angle", "[Utility.Tensor]") {
    const vec A{30, -20, 10, 30, -10, 20};

    const auto B = tensor::stress::to_tensor(tensor::dev(A));

    REQUIRE(tensor::stress::lode(eig_sym(B / sqrt(accu(square(B)))).eval()) == Approx(tensor::stress::lode(A)));
}

TEST_CASE("Norm of Stress", "[Utility.Tensor]") {
    const vec stress = {1., 2., 2., 2., 2., 0.};

    REQUIRE(Approx(5) == tensor::stress::norm(stress));
}

TEST_CASE("Norm of Strain", "[Utility.Tensor]") {
    const vec strain = {1., 2., 2., 4., 4., 0.};

    REQUIRE(Approx(5) == tensor::strain::norm(strain));
}

TEST_CASE("Rotation of Strain", "[Utility.Tensor]") {
    const auto strain = transform::strain::rotate({.01, -.01, 0.}, .25 * datum::pi);

    REQUIRE(Approx(.02) == tensor::strain::norm(strain));
}
