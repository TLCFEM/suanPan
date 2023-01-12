#include "CatchHeader.h"
#include <Toolbox/Expression.h>
#include <Toolbox/utility.h>

TEST_CASE("Variable Split", "[Utility.Expression]") {
    for(const auto test_list = std::vector<std::string>{
            "x|y|d_z",
            "\"x|y|d_z\"",
            R"("x""y""d_z")"
        }; const auto& I : test_list) {
        const auto variable_list = suanpan::expression::split(I);

        REQUIRE(variable_list.size() == 3);
        REQUIRE(variable_list[0] == "x");
        REQUIRE(variable_list[1] == "y");
        REQUIRE(variable_list[2] == "d_z");
    }
}

TEST_CASE("Expression Evaluation", "[Utility.Expression]") {
    auto expression = Expression(0, "x|y");

    REQUIRE(expression.compile("x^2+y^2+2*x*y") == true);

    mat test_data = randn(2, 100);

    suanpan_for(0llu, test_data.n_cols, [&](const uword I) {
        const auto& x = test_data(0, I);
        const auto& y = test_data(1, I);
        const auto f = x * x + y * y + 2. * x * y;
        REQUIRE(expression.evaluate(test_data.col(I)) == Approx(f));
        const auto df = 2. * (x + y);
        const auto gradient = expression.gradient(test_data.col(I));
        for(auto J = 0llu; J < test_data.n_rows; ++J)
            REQUIRE(gradient(J) == Approx(df));
    });
}
