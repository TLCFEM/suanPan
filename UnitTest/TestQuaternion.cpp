#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#include "CatchHeader.h"

#include <Toolbox/Quaternion.hpp>

TEST_CASE("Quaternion Basic Function", "[Utility.Quaternion]") {
    const Quaternion A(2., 3., 4., 5.);
    const Quaternion B(1., -2., 6., 3.);

    auto C = A + B;

    REQUIRE(Approx(3) == C.real());

    C = A * B;

    REQUIRE(Approx(-31) == C.real());
    REQUIRE(Approx(-19) == C.imag()(0));
    REQUIRE(Approx(-3) == C.imag()(1));
    REQUIRE(Approx(37) == C.imag()(2));

    C = B * A;

    REQUIRE(Approx(-31) == C.real());
    REQUIRE(Approx(17) == C.imag()(0));
    REQUIRE(Approx(35) == C.imag()(1));
    REQUIRE(Approx(-15) == C.imag()(2));

    C = B.inv() * B;

    const vec D(3, fill::randn);

    REQUIRE(norm(C * D - D) == Approx(0.));
}

TEST_CASE("Quaternion Conversion To Rotation Matrix", "[Utility.Quaternion]") {
    for(auto I = 0; I < 100; ++I) {
        Quaternion A(2., vec(4. * randu(3)));

        A.normalise();

        REQUIRE(A == transform::to_quaternion(A.to_mat()));
    }
}

TEST_CASE("Quaternion Conversion To Rotation Vector", "[Utility.Quaternion]") {
    for(auto I = 0; I < 100; ++I) {
        Quaternion A(as_scalar(randn(1)), vec(randu(3)));

        A.normalise();

        REQUIRE(A == transform::to_quaternion(A.to_pseudo()));
    }
}

TEST_CASE("Rodrigues Rotation", "[Utility.Rodrigues]") {
    for(auto I = 0; I < 100; ++I) {
        const auto A = transform::to_quaternion(vec(3, fill::randn));

        const auto B = transform::to_quaternion(transform::to_pseudo(A.to_mat()));

        REQUIRE(norm((B.inv() * A).to_pseudo()) < 1E-10);
    }
}

TEST_CASE("Quaternion Subtraction", "[Utility.Quaternion]") {
    const Quaternion A(2., 3., 4., 5.);
    const Quaternion B(1., 2., 3., 4.);

    const auto C = A - B;
    REQUIRE(Approx(1.) == C.real());
    REQUIRE(Approx(1.) == C.imag()(0));
    REQUIRE(Approx(1.) == C.imag()(1));
    REQUIRE(Approx(1.) == C.imag()(2));

    const auto D = A - A;
    REQUIRE(Approx(0.) == D.norm());
}

TEST_CASE("Quaternion Division", "[Utility.Quaternion]") {
    const Quaternion A(2., 3., 4., 5.);
    const Quaternion B(1., -2., 6., 3.);

    // A / B == A * inv(B)
    const auto C = A / B;
    const auto D = A * B.inv();
    REQUIRE(Approx(C.real()).margin(1E-12) == D.real());
    REQUIRE(norm(C.imag() - D.imag()) == Approx(0.).margin(1E-12));

    // |A / B| == |A| / |B|
    REQUIRE(Approx(std::sqrt(C.norm())) == std::sqrt(A.norm()) / std::sqrt(B.norm()));

    const auto E = B / B;
    REQUIRE(Approx(1.) == E.real());
    REQUIRE(norm(E.imag()) == Approx(0.).margin(1E-12));
}

TEST_CASE("Quaternion Conjugate And Norm", "[Utility.Quaternion]") {
    const Quaternion A(1., 2., 3., 4.);

    const auto c = A.conj();
    REQUIRE(Approx(1.) == c.real());
    REQUIRE(Approx(-2.) == c.imag()(0));
    REQUIRE(Approx(-3.) == c.imag()(1));
    REQUIRE(Approx(-4.) == c.imag()(2));

    const auto h = A * A.conj();
    REQUIRE(Approx(1. + 4. + 9. + 16.) == h.real());
    REQUIRE(norm(h.imag()) == Approx(0.).margin(1E-12));
    REQUIRE(Approx(1. + 4. + 9. + 16.) == A.norm());
}

TEST_CASE("Quaternion Compound Assignment", "[Utility.Quaternion]") {
    Quaternion A(1., 2., 3., 4.);
    const Quaternion B(5., 6., 7., 8.);

    A += B;
    REQUIRE(Approx(6.) == A.real());
    A -= B;
    REQUIRE(Approx(1.) == A.real());

    A *= B;
    const auto expected = Quaternion(1., 2., 3., 4.) * B;
    REQUIRE(Approx(expected.real()) == A.real());
    REQUIRE(Approx(expected.imag()(0)) == A.imag()(0));

    A /= B;
    REQUIRE(Approx(1.).margin(1E-12) == A.real());
    REQUIRE(Approx(2.).margin(1E-12) == A.imag()(0));
}

TEST_CASE("Quaternion Unary Minus And Rotating Vector", "[Utility.Quaternion]") {
    const Quaternion A(1., 2., 3., 4.);
    const auto neg = -A;
    REQUIRE(Approx(-1.) == neg.real());

    // rotating a vector by unit quaternion preserves its length
    Quaternion Q(.5, vec{1., 0., 0.});
    Q.normalise();
    const vec v = {1., 2., 3.};
    const auto rotated = Q * v;
    REQUIRE(norm(rotated) == Approx(norm(v)).margin(1E-12));
}
#pragma GCC diagnostic pop
