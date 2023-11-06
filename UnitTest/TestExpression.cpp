#include "CatchHeader.h"
#include <Toolbox/Expression.h>
#include <Toolbox/utility.h>

TEST_CASE("Variable Split", "[Utility.Expression]") {
    for(const auto test_list = std::vector<std::string>{"x|y|d_z", "\"x|y|d_z\"", R"("x""y""d_z")"}; const auto& I : test_list) {
        const auto variable_list = suanpan::expression::split(I);

        REQUIRE(variable_list.size() == 3);
        REQUIRE(variable_list[0].first == "x");
        REQUIRE(variable_list[1].first == "y");
        REQUIRE(variable_list[2].first == "d_z");
    }
}

TEST_CASE("Simple Scalar Evaluation", "[Utility.Expression]") {
    const auto expression = make_unique<SimpleScalarExpression>(0, "x|y");

    REQUIRE(expression->compile("x^2+y^2+2*x*y") == true);

    mat test_data = randu(2, 100);

    suanpan_for(0llu, test_data.n_cols, [&](const uword I) {
        const auto expression_copy = expression->get_copy();
        const auto &x = test_data(0, I), &y = test_data(1, I);
        const auto f = x * x + y * y + 2. * x * y;
        REQUIRE(expression_copy->evaluate(test_data.col(I)).at(0) == Approx(f));
        const auto df = 2. * (x + y);
        const auto gradient = expression_copy->gradient(test_data.col(I));
        for(auto J = 0llu; J < test_data.n_rows; ++J)
            REQUIRE(gradient(J) == Approx(df));
    });
}

TEST_CASE("Simple Dot Evaluation", "[Utility.Expression]") {
    const auto expression = make_unique<SimpleScalarExpression>(0, "x|2");

    REQUIRE(expression->compile("sum(x)") == true);

    mat test_data = randu(2, 100);

    suanpan_for(0llu, test_data.n_cols, [&](const uword I) {
        const auto expression_copy = expression->get_copy();
        const auto &x = test_data(0, I), &y = test_data(1, I);
        const auto f = x + y;
        REQUIRE(expression_copy->evaluate(test_data.col(I)).at(0) == Approx(f));
        const auto gradient = expression_copy->gradient(test_data.col(I));
        for(auto J = 0llu; J < test_data.n_rows; ++J)
            REQUIRE(gradient(J) == Approx(1));
    });
}

TEST_CASE("Simple Vector Evaluation", "[Utility.Expression]") {
    const auto expression = make_unique<SimpleVectorExpression>(0, "x|3", "\"y|2\"");

    REQUIRE(expression->compile("y[0]:=x[0]+x[1]+x[2];y[1]:=x[0]*x[1]*x[2];") == true);

    mat test_data = randu(3, 100);

    suanpan_for(0llu, test_data.n_cols, [&](const uword I) {
        const auto expression_copy = expression->get_copy();
        const auto f = expression_copy->evaluate(test_data.col(I));
        REQUIRE(f.at(0) == Approx(sum(test_data.col(I))));
        REQUIRE(f.at(1) == Approx(prod(test_data.col(I))));
    });
}
