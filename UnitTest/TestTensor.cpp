#include <Toolbox/tensor.h>
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

TEST_CASE("Basic Quantities", "[Utility.Tensor]") {
    const vec A = dev(dev(unit_symmetric_tensor4())).diag();
    const auto B = stress::norm(A), C = stress::norm(vec(A)), D = strain::norm(A), E = strain::norm(vec(A));

    REQUIRE(Approx(B) == C);
    REQUIRE(Approx(D) == E);

    REQUIRE(norm(strain::to_green(eye(3, 3))) == Approx(0));

    REQUIRE(norm(transform::compute_jacobian_nominal_to_principal(strain::to_green(eye(2, 2)))) == Approx(0));

    for(auto I = 0; I < 1000; ++I) {
        vec F(6, fill::randn);

        REQUIRE(norm(strain::to_voigt(strain::to_tensor(F)) - F) <= 1E-13);
        REQUIRE(norm(stress::to_voigt(stress::to_tensor(F)) - F) <= 1E-13);

        F = randn(3);

        REQUIRE(norm(strain::to_voigt(strain::to_tensor(F)) - F) <= 1E-13);
        REQUIRE(norm(stress::to_voigt(stress::to_tensor(F)) - F) <= 1E-13);
    }

    transform::compute_jacobian_principal_to_nominal(randn(2, 2));
}

TEST_CASE("Base Conversion", "[Utility.Tensor]") {
    const vec3 g1{2, 0, 0}, g2{1, 1, 0}, g3{1, 2, 3};
    const auto g = base::Base3D(g1, g2, g3);
    const auto [g1_, g2_, g3_] = g.to_inverse();

    REQUIRE(Approx(g1_(0)).margin(1E-15) == .5);
    REQUIRE(Approx(g1_(1)).margin(1E-15) == -.5);
    REQUIRE(Approx(g1_(2)).margin(1E-15) == 1. / 6.);
    REQUIRE(Approx(g2_(0)).margin(1E-15) == 0.);
    REQUIRE(Approx(g2_(1)).margin(1E-15) == 1.);
    REQUIRE(Approx(g2_(2)).margin(1E-15) == -2. / 3.);
    REQUIRE(Approx(g3_(0)).margin(1E-15) == 0.);
    REQUIRE(Approx(g3_(1)).margin(1E-15) == 0.);
    REQUIRE(Approx(g3_(2)).margin(1E-15) == 1. / 3.);
}

TEST_CASE("Curvature Tensor", "[Utility.Tensor]") {
    for(auto I = 0; I < 100; ++I) {
        const auto z = randn(2);
        const vec3 a1{1., 0., z(0)}, a2{0., 1., -z(1)};
        const auto a3 = base::unit_norm(a1, a2);
        const auto t = 1. / sqrt(1. + accu(square(z)));
        mat22 b(fill::zeros);
        b(0, 0) = dot(a3, vec{0., 0., 1.});
        b(1, 1) = dot(a3, vec{0., 0., -1.});
        REQUIRE(Approx(b(0, 0)).margin(1E-15) == t);
        REQUIRE(Approx(b(1, 1)).margin(1E-15) == -t);
        mat22 a;
        a(0, 0) = dot(a1, a1);
        a(1, 1) = dot(a2, a2);
        a(0, 1) = dot(a1, a2);
        a(1, 0) = a(0, 1);
        const mat22 c = b * solve(a, b);
        REQUIRE(Approx(c(0, 0)).margin(1E-15) == (1. + z(1) * z(1)) * pow(t, 4.));
        REQUIRE(Approx(c(1, 1)).margin(1E-15) == (1. + z(0) * z(0)) * pow(t, 4.));
        REQUIRE(Approx(c(0, 1)).margin(1E-15) == -z(0) * z(1) * pow(t, 4.));
        REQUIRE(Approx(c(1, 0)).margin(1E-15) == -z(0) * z(1) * pow(t, 4.));
    }
}
