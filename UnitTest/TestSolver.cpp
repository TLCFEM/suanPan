#include "TestSolver.h"
#include <Toolbox/LBFGS.hpp>
#include "CatchHeader.h"

TEST_CASE("LBFGS Solver", "[Utility.Solver]") {
    Quadratic function;

    vec result{2., 3.};

    LBFGS().optimize(function, result);

    for(const auto I : result)
        REQUIRE(Approx(1.) == I);
}
