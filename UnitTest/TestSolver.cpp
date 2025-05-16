#include "CatchHeader.h"

#include <Toolbox/LBFGS.hpp>
#include <Toolbox/brent.hpp>
#include <Toolbox/ridders.hpp>

class Quadratic {
public:
    [[nodiscard]] static vec evaluate_residual(const vec& x) { return square(x) - 1.; }

    [[nodiscard]] static mat evaluate_jacobian(const vec& x) { return diagmat(2. * x); }
};

TEST_CASE("LBFGS Solver", "[Utility.Solver]") {
    Quadratic function;

    vec result{2., 3.};

    LBFGS().optimize(function, result);

    for(const auto I : result) REQUIRE_THAT(I, Catch::Matchers::WithinAbs(1., 100 * datum::eps));
}

TEST_CASE("Ridders Solver", "[Utility.Solver]") {
    auto fn = [](const double x) { return std::pow(x - 1., 3); };

    const auto root = ridders(fn, 0., fn(0.), 4., fn(4.), datum::eps);

    REQUIRE_THAT(fn(root), Catch::Matchers::WithinAbs(0., datum::eps));
}

TEST_CASE("Brent Solver", "[Utility.Solver]") {
    auto fn = [](const double x) { return std::pow(x - 1., 3); };

    const auto root = brent(fn, 0., 4., 1E-12);

    REQUIRE_THAT(fn(root), Catch::Matchers::WithinAbs(0., 1E-12));
}
