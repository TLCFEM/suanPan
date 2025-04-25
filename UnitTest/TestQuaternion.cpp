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
#pragma GCC diagnostic pop
