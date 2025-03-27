#include "CatchHeader.h"

#include <Toolbox/LBFGS.hpp>

class Quadratic {
public:
    [[nodiscard]] static vec evaluate_residual(const vec& x) { return square(x) - 1.; }

    [[nodiscard]] static mat evaluate_jacobian(const vec& x) { return diagmat(2. * x); }
};

TEST_CASE("LBFGS Solver", "[Utility.Solver]") {
    Quadratic function;

    vec result{2., 3.};

    LBFGS().optimize(function, result);

    for(const auto I : result)
        REQUIRE(Approx(1.) == I);
}
