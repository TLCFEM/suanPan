#include "CatchHeader.h"

#include <Toolbox/homotopy.hpp>

TEST_CASE("Arctan", "[Utility.Homotopy]") {
    const auto arctan_system = [](const vec& x) {
        return std::pair{vec{std::atan(x(0))}, mat{1. / (1. + x(0) * x(0))}};
    };

    vec x{3.};

    homotopy_config config;
    config.incre_t = .1;
    config.tolerance = 1e-10;
    config.min_incre = 1e-3;
    config.max_evaluation = 200u;
    config.max_iteration = 10u;

    const int status = homotopy_solve<mat>(x, arctan_system, config);

    REQUIRE(status == SUANPAN_SUCCESS);
    REQUIRE_THAT(x(0), Catch::Matchers::WithinAbs(0., 1e-9));
}
