#include "CatchHeader.h"

#include <Toolbox/tensor.h>

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

TEST_CASE("Suanpan Ramp", "[Utility.Tensor]") {
    REQUIRE(Approx(0.) == suanpan::ramp(-3.));
    REQUIRE(Approx(0.) == suanpan::ramp(0.));
    REQUIRE(Approx(2.5) == suanpan::ramp(2.5));
}

TEST_CASE("Tensor Trace", "[Utility.Tensor]") {
    const vec a{1., 2., 3., 4., 5., 6.};
    REQUIRE(Approx(6.) == tensor::trace<3>(a));
    REQUIRE(Approx(-6.) == tensor::trace<3>(-a));
    REQUIRE(Approx(21.) == tensor::trace<6>(a));
    REQUIRE(Approx(-21.) == tensor::trace<6>(-a));
    // mean is trace divided by dimension
    REQUIRE(Approx(2.) == tensor::mean<3>(a));
    REQUIRE(Approx(3.5) == tensor::mean<6>(a));
}

TEST_CASE("Tensor Deviatoric", "[Utility.Tensor]") {
    const vec a{3., 6., 9., 2., 4., 8.};
    const auto m = tensor::mean<3>(a);
    const auto d = dev(a);
    REQUIRE(Approx(0.) == d(0) - (3. - m));
    REQUIRE(Approx(0.) == d(1) - (6. - m));
    REQUIRE(Approx(0.) == d(2) - (9. - m));
    // first three deviatoric components must sum to zero
    REQUIRE(Approx(0.) == sum(d.head(3)));
}

TEST_CASE("Tensor Unit Tensors", "[Utility.Tensor]") {
    const auto I = unit_symmetric_tensor4();
    REQUIRE(I.n_rows == I.n_cols);
    REQUIRE(6 == I.n_rows);
    // the symmetric fourth-order identity keeps normal (Voigt 0..2) components
    // and halves shear (Voigt 3..5) components (Voigt factor of 2)
    const vec strain = {1., 2., 3., 4., 5., 6.};
    const vec out = I * strain;
    for(auto J = 0u; J < 3u; ++J) REQUIRE(Approx(strain(J)) == out(J));
    for(auto J = 3u; J < 6u; ++J) REQUIRE(Approx(.5 * strain(J)) == out(J));

    const auto id = unit_deviatoric_tensor4();
    // projecting a purely volumetric strain (equal normal, zero shear in Voigt) yields zero deviator
    const vec volumetric = {1., 1., 1., 0., 0., 0.};
    const auto devproj = id * volumetric;
    REQUIRE(Approx(0.).margin(1E-12) == accu(devproj));

    // projecting a purely deviatoric normal strain (trace-free, no shear) reproduces it
    const vec deviatoric = {1., 1., -2., 0., 0., 0.};
    const auto devproj2 = id * deviatoric;
    REQUIRE(norm(devproj2 - deviatoric) == Approx(0.).margin(1E-12));

    const auto id2 = unit_deviatoric_tensor4v2();
    REQUIRE(id2.n_rows == id2.n_cols);
    REQUIRE(6 == id2.n_rows);
}

TEST_CASE("Tensor Isotropic Stiffness", "[Utility.Tensor]") {
    constexpr double E = 200E3, NU = .3;
    const auto D = isotropic_stiffness(E, NU);
    REQUIRE(6 == D.n_rows);
    REQUIRE(6 == D.n_cols);
    // symmetric modulus tensor
    REQUIRE(Approx(0.).margin(1E-6) == accu(D - D.t()));

    constexpr double G = E / (2. + 2. * NU);
    constexpr double lambda = G * NU / (.5 - NU);
    // verify a few closed-form entries
    REQUIRE(Approx(2. * G + lambda) == D(0, 0));
    REQUIRE(Approx(lambda) == D(0, 1));
    REQUIRE(Approx(G) == D(3, 3));

    // extra dimension padding
    const auto D2 = isotropic_stiffness(E, NU, 2);
    REQUIRE(8 == D2.n_rows);
    REQUIRE(8 == D2.n_cols);
    const auto D3 = isotropic_stiffness(E, NU, 3);
    REQUIRE(9 == D3.n_rows);
    REQUIRE(9 == D3.n_cols);
}

TEST_CASE("Tensor Orthotropic Stiffness", "[Utility.Tensor]") {
    const vec mod = {100., 50., 30., 20., 10., 15.};
    const vec ratio = {.2, .3, .1};
    const auto D = orthotropic_stiffness(mod, ratio);
    REQUIRE(6 == D.n_rows);
    REQUIRE(6 == D.n_cols);
    // symmetric compliance inverse construction => symmetric modulus
    REQUIRE(Approx(0.).margin(1E-8) == accu(D - D.t()));
    // shear moduli appear verbatim on the diagonal
    REQUIRE(Approx(20.) == D(3, 3));
    REQUIRE(Approx(10.) == D(4, 4));
    REQUIRE(Approx(15.) == D(5, 5));
}

TEST_CASE("Skew Symmetric Matrix", "[Utility.Tensor]") {
    const vec r{1., 2., 3.};
    const auto S = transform::skew_symm(r);
    REQUIRE(3 == S.n_rows);
    REQUIRE(3 == S.n_cols);
    // skew symmetric => S = -S.t()
    REQUIRE(Approx(0.) == accu(S + S.t()));
    REQUIRE(Approx(0.) == S(0, 0));
    REQUIRE(Approx(0.) == S(1, 1));
    REQUIRE(Approx(0.) == S(2, 2));
    REQUIRE(Approx(-3.) == S(0, 1));
    REQUIRE(Approx(3.) == S(1, 0));
    REQUIRE(Approx(-1.) == S(1, 2));
    REQUIRE(Approx(1.) == S(2, 1));
    REQUIRE(Approx(2.) == S(0, 2));
    REQUIRE(Approx(-2.) == S(2, 0));
}

TEST_CASE("Rodrigues Rotation", "[Utility.Tensor]") {
    // a rotation of pi about z axis maps (1,0) to (-1,0)
    const vec axis{0., 0., datum::pi};
    const auto R = transform::rodrigues(axis);
    REQUIRE(3 == R.n_rows);
    REQUIRE(3 == R.n_cols);
    REQUIRE(norm(R * vec{1., 0., 0.} - vec{-1., 0., 0.}) == Approx(0.).margin(1E-9));

    // identity rotation about the zero axis
    const vec zero{0., 0., 0.};
    REQUIRE(norm(transform::rodrigues(zero) - eye(3, 3)) == Approx(0.).margin(1E-12));

    // rotation matrices must be orthogonal
    for(auto I = 0; I < 100; ++I) {
        const auto RR = transform::rodrigues(randn(3));
        REQUIRE(norm(RR * RR.t() - eye(3, 3)) == Approx(0.).margin(1E-10));
        REQUIRE(Approx(1.).margin(1E-9) == det(RR));
    }
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
