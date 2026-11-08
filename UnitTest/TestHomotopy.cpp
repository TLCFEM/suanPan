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

TEST_CASE("Quadratic System", "[Utility.Homotopy]") {
    // solve x^2 - 4 = 0 returning Jacobian 2x
    const auto system = [](const vec& x) {
        return std::pair{vec{x(0) * x(0) - 4.}, mat{2. * x(0)}};
    };

    vec x{1.}; // must converge to the positive root 2

    homotopy_config config;
    config.incre_t = .1;
    config.tolerance = 1e-10;

    const int status = homotopy_solve<mat>(x, system, config);

    REQUIRE(status == SUANPAN_SUCCESS);
    REQUIRE_THAT(x(0), Catch::Matchers::WithinAbs(2., 1e-8));
}

TEST_CASE("Linear System Multi Variable", "[Utility.Homotopy]") {
    // 2-variable decoupled linear system:
    //   2 x0 = 2
    //   2 x1 = 4  => solution x0=1, x1=2
    const auto system = [](const vec& x) {
        mat J(2, 2, fill::zeros);
        J(0, 0) = 2.;
        J(1, 1) = 2.;
        vec f(2, fill::none);
        f(0) = 2. * x(0) - 2.;
        f(1) = 2. * x(1) - 4.;
        return std::pair{std::move(f), std::move(J)};
    };

    vec x{0., 0.};

    homotopy_config config;
    config.incre_t = .2;
    config.tolerance = 1e-10;
    config.fixed_incre = true;

    const int status = homotopy_solve<mat>(x, system, config);

    REQUIRE(status == SUANPAN_SUCCESS);
    REQUIRE_THAT(x(0), Catch::Matchers::WithinAbs(1., 1e-8));
    REQUIRE_THAT(x(1), Catch::Matchers::WithinAbs(2., 1e-8));
}
