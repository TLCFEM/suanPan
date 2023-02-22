#include <Toolbox/IntegrationPlan.h>
#include "CatchHeader.h"

double test_integration(const IntegrationPlan& plan, double (*func)(double)) {
    double result = 0.;
    for(auto J = 0u; J < plan.n_rows; ++J) result += func(plan(J, 0)) * plan(J, 1);
    return result;
}

double test_cubic(const double x) { return x * x * x + 2.; }

double test_quartic(const double x) { return x * x * x * x - 2. * x * x * x - 6. * x * x + 2. * x - 1.; }

TEST_CASE("Radau", "[Utility.Integration]") {
    for(auto I = 3u; I <= 10u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::RADAU), test_cubic)) == 4);
    for(auto I = 4u; I <= 10u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::RADAU), test_quartic)) == -5.6);
}

TEST_CASE("Lobatto", "[Utility.Integration]") {
    for(auto I = 3u; I <= 20u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::LOBATTO), test_cubic)) == 4);
    for(auto I = 4u; I <= 20u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::LOBATTO), test_quartic)) == -5.6);
}

TEST_CASE("Gauss", "[Utility.Integration]") {
    for(auto I = 2u; I <= 20u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::GAUSS), test_cubic)) == 4);
    for(auto I = 3u; I <= 20u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::GAUSS), test_quartic)) == -5.6);
}

TEST_CASE("Chebyshev", "[Utility.Integration]") {
    for(auto I = 2u; I <= 7u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::CHEBYSHEV), test_cubic)) == 4);
    for(auto I = 4u; I <= 7u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::CHEBYSHEV), test_quartic)) == -5.6);
}

double test_weighted_function(const double x) { return x * x + 2.; }

TEST_CASE("Hermite", "[Utility.Integration]") {
    for(auto I = 2u; I <= 20u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::HERMITE), test_weighted_function)) == 4.43113);
}

TEST_CASE("Laguerre", "[Utility.Integration]") {
    for(auto I = 2u; I <= 20u; ++I)
        REQUIRE(Approx(test_integration(IntegrationPlan(1, I, IntegrationType::LAGUERRE), test_weighted_function)) == 4);
}
