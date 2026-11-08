#include "CatchHeader.h"

#include <Toolbox/brent.hpp>
#include <Toolbox/ridders.hpp>
#include <cmath>

TEST_CASE("Brent Linear Root", "[Utility.RootFinding]") {
    // f(x) = 2x - 6 has root x=3
    const auto f = [](const double x) { return 2. * x - 6.; };
    const auto root = brent(f, 0., 10., 1E-10);
    REQUIRE_THAT(root, Catch::Matchers::WithinAbs(3., 1E-8));
}

TEST_CASE("Brent Quadratic Root", "[Utility.RootFinding]") {
    // f(x) = x^2 - 4 has roots x=+-2
    const auto f = [](const double x) { return x * x - 4.; };
    const auto root = brent(f, 0., 5., 1E-10);
    REQUIRE_THAT(root, Catch::Matchers::WithinAbs(2., 1E-8));
}

TEST_CASE("Brent Sinusoidal Root", "[Utility.RootFinding]") {
    // sin(x) has a root at pi within [3, 4]
    const auto f = [](const double x) { return std::sin(x); };
    const auto root = brent(f, 3., 4., 1E-12);
    REQUIRE_THAT(root, Catch::Matchers::WithinAbs(datum::pi, 1E-8));
}

TEST_CASE("Brent Float Type", "[Utility.RootFinding]") {
    // check the float template instantiation works
    const auto f = [](const float x) { return x - 2.5f; };
    const auto root = brent(f, 0.f, 10.f, 1E-5f);
    REQUIRE_THAT(root, Catch::Matchers::WithinAbs(2.5f, 1E-3f));
}

TEST_CASE("Ridders Quadratic Root", "[Utility.RootFinding]") {
    // f(x) = x^2 - 9 has root x=3
    const auto f = [](const double x) { return x * x - 9.; };
    const auto root = ridders(f, 0., 10., 1E-10);
    REQUIRE_THAT(root, Catch::Matchers::WithinAbs(3., 1E-8));
}

TEST_CASE("Ridders Cos Root", "[Utility.RootFinding]") {
    // cos(x) has a root at pi/2 within [0, 2]
    const auto f = [](const double x) { return std::cos(x); };
    const auto root = ridders(f, 0., 2., 1E-12);
    REQUIRE_THAT(root, Catch::Matchers::WithinAbs(.5 * datum::pi, 1E-8));
}

TEST_CASE("Ridders Guess", "[Utility.RootFinding]") {
    // ridders_guess expands the interval until a sign change is found
    const auto f = [](const double x) { return x - 4.; };
    // start with an interval that does not bracket the root
    const auto root = ridders_guess(f, 0., 1., 1E-10);
    REQUIRE_THAT(root, Catch::Matchers::WithinAbs(4., 1E-8));
}

TEST_CASE("Ridders Guess With Data", "[Utility.RootFinding]") {
    // exp(x) - 100 has root x = ln(100) ~ 4.605
    const auto f = [](const double x) { return std::exp(x) - 100.; };
    constexpr double x1 = 0., x2 = 1.;
    const auto root = ridders_guess(f, x1, f(x1), x2, f(x2), 1E-10);
    REQUIRE_THAT(root, Catch::Matchers::WithinAbs(std::log(100.), 1E-8));
}
