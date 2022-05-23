#include <Toolbox/tensorToolbox.h>
#include "CatchHeader.h"

using namespace tensor;

TEST_CASE("Fixed Stress Invariant", "[Utility.Tensor]") {
    const vec A{50, -20, 10, 30, -10, 20};

    REQUIRE(Approx(40) == tensor::stress::invariant1(A));
    REQUIRE(Approx(2100) == tensor::stress::invariant2(A));
    REQUIRE(Approx(-28000) == tensor::stress::invariant3(A));

    const auto B = eig_sym(stress::to_tensor(A)).eval();

    REQUIRE(Approx(40) == tensor::stress::invariant1(B));
    REQUIRE(Approx(2100) == tensor::stress::invariant2(B));
    REQUIRE(Approx(-28000) == tensor::stress::invariant3(B));
}

TEST_CASE("Random Stress/Strain Invariant", "[Utility.Tensor]") {
    for(auto I = 0; I < 10000; ++I) {
        const vec A(6, fill::randn);

        mat M = stress::to_tensor(A), E;
        vec B;

        eig_sym(B, E, M, "std");

        REQUIRE(std::fabs(stress::invariant1(A) - stress::invariant1(B)) <= 1E-13);
        REQUIRE(std::fabs(stress::invariant2(A) - stress::invariant2(B)) <= 1E-13);
        REQUIRE(std::fabs(stress::invariant3(A) - stress::invariant3(B)) <= 1E-13);

        M = strain::to_tensor(A);

        eig_sym(B, E, M, "std");

        REQUIRE(std::fabs(strain::invariant1(A) - strain::invariant1(B)) <= 1E-13);
        REQUIRE(std::fabs(strain::invariant2(A) - strain::invariant2(B)) <= 1E-13);
        REQUIRE(std::fabs(strain::invariant3(A) - strain::invariant3(B)) <= 1E-13);
    }
}

TEST_CASE("Lode Angle", "[Utility.Tensor]") {
    const vec A{30, -20, 10, 30, -10, 20};

    const auto B = stress::to_tensor(dev(A));

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

TEST_CASE("Basic Quantites", "[Utility.Tensor]") {
    const vec A = dev(dev(unit_symmetric_tensor4())).diag();
    const auto B = stress::norm(A), C = stress::norm(vec(A)), D = strain::norm(A), E = strain::norm(vec(A));

    REQUIRE(Approx(B) == C);
    REQUIRE(Approx(D) == D);
}
