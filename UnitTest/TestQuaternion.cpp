#include <catch/catch.hpp>
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
