#include <Element/Utility/IGA/NURBSSurface.h>
#include <Element/Utility/IGA/NURBSVolume.h>
#include "CatchHeader.h"

TEST_CASE("BSpline Compute Span", "[IGA.BSpline]") {
    const auto A = BSplineCurve2D(vec{0., 0., 0., .1, .3, .7, 1., 1., 1.});

    REQUIRE(5 == A.evaluate_span(.9));
    REQUIRE(2 == A.evaluate_span(.0));
    REQUIRE(5 == A.evaluate_span(1.));
}

TEST_CASE("BSpline Compute Basis Function", "[IGA.BSpline]") {
    const auto A = BSplineCurve2D(vec{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5});

    std::random_device rd;
    std::mt19937 gen(rd());
    // ReSharper disable once CppLocalVariableMayBeConst
    std::uniform_real_distribution dis(2.0, 3.0);

    const auto B = dis(gen);

    mat C = A.evaluate_basis(B);

    REQUIRE(C(0) == Approx(.5 * (3. - B) * (3. - B)));
    REQUIRE(C(1) == Approx(-5.5 + (5. - B) * B));
    REQUIRE(C(2) == Approx(.5 * (B - 2.) * (B - 2.)));

    C = A.evaluate_basis_derivative(B, 4);

    REQUIRE(C(1, 0) == Approx(B - 3.));
    REQUIRE(C(1, 1) == Approx(5. - 2. * B));
    REQUIRE(C(1, 2) == Approx(B - 2.));
    REQUIRE(C(2, 0) == Approx(1.));
    REQUIRE(C(2, 1) == Approx(-2.));
    REQUIRE(C(2, 2) == Approx(1.));
}

TEST_CASE("BSpline Compute Derivative of Point", "[IGA.BSpline]") {
    const auto A = BSplineCurve3D(vec{0, 0, 0, 1, 2, 3, 3, 3});

    field<vec> P(5);
    P.for_each([](vec& t_point) { t_point.randu(3); });

    SECTION("Point Derivative Matches Shape Function") {
        const auto C = A.evaluate_point_derivative(2.5, P, 2)(2);
        const auto DS = A.evaluate_shape_function_derivative(2.5, 2)(2);

        vec D(3, fill::zeros);

        for(auto J = 0; J < 5; ++J) D += P(J) * DS(J);

        REQUIRE(norm(C) == Approx(norm(D)).epsilon(1E-10));
    }
}

TEST_CASE("BSplineSurface Compute Location of A Point", "[IGA.BSplineSurface]") {
    const auto U = BSplineSurface4D(vec{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}, vec{0, 0, 0, 1, 2, 3, 3, 3});

    field<vec> P(8, 5);
    P(2, 1) = vec{0, 2, 4, 1};
    P(2, 2) = vec{0, 6, 4, 2};
    P(2, 3) = vec{0, 2, 0, 1};
    P(3, 1) = vec{4, 6, 8, 2};
    P(3, 2) = vec{12, 24, 12, 6};
    P(3, 3) = vec{4, 6, 0, 2};
    P(4, 1) = vec{4, 2, 4, 1};
    P(4, 2) = vec{8, 6, 4, 2};
    P(4, 3) = vec{4, 2, 0, 1};

    const auto A = U.evaluate_point(2.5, 1., P);

    REQUIRE(A(0) == Approx(54. / 8.));
    REQUIRE(A(1) == Approx(98. / 8.));
    REQUIRE(A(2) == Approx(68. / 8.));
    REQUIRE(A(3) == Approx(27. / 8.));

    SECTION("Matches Shape Function") {
        const auto B = U.evaluate_shape_function(2.5, 1., P);

        vec PP = zeros(4);
        for(auto I = 0llu; I < B.n_rows; ++I) for(auto J = 0llu; J < B.n_cols; ++J) if(!P(I, J).empty()) PP += B(I, J) * P(I, J);

        REQUIRE(norm(PP) == Approx(norm(A)));
    }
}

TEST_CASE("BSplineSurface Compute Derivative of A Point", "[IGA.BSplineSurface]") {
    const auto U = BSplineSurface4D(vec{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}, vec{0, 0, 0, 1, 2, 3, 3, 3});

    field<vec> P(10, 10);
    P.for_each([](vec& N) { N.randu(4); });

    const auto A = U.evaluate_point_derivative(0., 0., P, 3, 3);

    SECTION("Compute Derivative of A Point") {
        const vec B = A(1, 0) - 2. * (P(1, 0) - P(0, 0));
        const vec C = A(0, 1) - 2. * (P(0, 1) - P(0, 0));
        const vec D = A(1, 1) - 4. * (P(1, 1) + P(0, 0) - P(0, 1) - P(1, 0));

        REQUIRE(norm(B) == Approx(0.));
        REQUIRE(norm(C) == Approx(0.));
    }

    SECTION("Matches Shape Function") {
        const auto E = U.evaluate_shape_function_derivative(0., 0., 1, 1)(1, 1);

        vec PP = zeros(4);
        for(auto I = 0llu; I < E.n_rows; ++I) for(auto J = 0llu; J < E.n_cols; ++J) PP += E(I, J) * P(I, J);

        REQUIRE(norm(A(1, 1)) == Approx(norm(PP)));
    }
}

TEST_CASE("BSplineVolume Compute Location of A Point", "[IGA.BSplineVolume]") {
    auto U = BSplineVolume4D(vec{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}, vec{0, 0, 0, 1, 2, 3, 3, 3}, vec{0, 0, 0, 1, 2, 3, 3, 3});

    field<vec> P(10, 10, 10);
    P.for_each([](vec& N) { N.randu(4); });

    U.set_control_polygon(P);

    const auto A = U.evaluate_point(2.5, 1., .5, P);

    const auto B = U.evaluate_shape_function(2.5, 1., .5);

    vec C = zeros(4);

    for(auto I = 0llu; I < B.n_rows; ++I) for(auto J = 0llu; J < B.n_cols; ++J) for(auto K = 0llu; K < B.n_slices; ++K) C += B(I, J, K) * P(I, J, K);

    REQUIRE(norm(A) == Approx(norm(C)));

    SECTION("Matches Shape Function") {
        const auto D = U.evaluate_point_derivative(2.5, 1., .5, 1)(1, 1, 1);

        const auto E = U.evaluate_shape_function_derivative(2.5, 1., .5, 1, 1, 1)(1, 1, 1);

        C.zeros(4);
        for(auto I = 0llu; I < E.n_rows; ++I) for(auto J = 0llu; J < E.n_cols; ++J) for(auto K = 0llu; K < E.n_slices; ++K) C += E(I, J, K) * P(I, J, K);

        REQUIRE(norm(D) == Approx(norm(C)));
    }
}

TEST_CASE("NURBS Compute Location of A Point", "[IGA.NURBS]") {
    const auto A = NURBSCurve2D(vec{0, 0, 0, 1, 2, 3, 3, 3});

    field<vec> P(5, 1);

    P(0) = {0, 0, 1};
    P(1) = {1, 1, 4};
    P(2) = {3, 2, 1};
    P(3) = {4, 1, 1};
    P(4) = {5, -1, 1};

    IGA::convert_to_weighted(P);

    const auto B = A.evaluate_point(1, P);

    REQUIRE(B(0) == Approx(1.4));
    REQUIRE(B(1) == Approx(1.2));

    const auto C = A.evaluate_shape_function(1., P);

    vec BB(3, fill::zeros);
    for(auto I = 0llu; I < P.n_rows; ++I) BB += P(I) * C(I);

    REQUIRE(B(0) == Approx(BB(0)));
    REQUIRE(B(1) == Approx(BB(1)));
}

// EX4.2
TEST_CASE("NURBS Compute Derivative of A Point", "[IGA.NURBS]") {
    const auto A = NURBSCurve2D(vec{0, 0, 0, 1, 1, 1});

    field<vec> P(3, 1);

    P(0) = {1, 0, 1};
    P(1) = {1, 1, 1};
    P(2) = {0, 1, 2};

    auto PP = P;

    IGA::convert_to_weighted(P);

    auto B = A.evaluate_point_derivative(0, P, 2);
    const auto C = A.evaluate_point_derivative(1, P, 6);

    REQUIRE(B(1)(0) == Approx(0.));
    REQUIRE(B(1)(1) == Approx(2.));
    REQUIRE(B(2)(0) == Approx(-4.));
    REQUIRE(B(2)(1) == Approx(0.));
    REQUIRE(C(1)(0) == Approx(-1.));
    REQUIRE(C(1)(1) == Approx(0.));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution dis(.0, 1.);

    const auto E = dis(gen);

    B = A.evaluate_point_derivative(E, P, 6);

    const auto D = A.evaluate_shape_function_derivative(E, PP, 4)(4);

    vec BB(3, fill::zeros);
    for(auto I = 0llu; I < PP.n_rows; ++I) if(!PP(I).empty()) BB += PP(I) * D(I);

    REQUIRE(norm(BB.head(2)) == Approx(norm(B(4))));
}

TEST_CASE("NURBSSurface Compute Location of A Point", "[IGA.NURBSSurface]") {
    const auto U = NURBSSurface3D(vec{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}, vec{0, 0, 0, 1, 2, 3, 3, 3});

    field<vec> P(5, 4);
    P(2, 1) = vec{0, 2, 4, 1};
    P(2, 2) = vec{0, 6, 4, 2};
    P(2, 3) = vec{0, 2, 0, 1};
    P(3, 1) = vec{4, 6, 8, 2};
    P(3, 2) = vec{12, 24, 12, 6};
    P(3, 3) = vec{4, 6, 0, 2};
    P(4, 1) = vec{4, 2, 4, 1};
    P(4, 2) = vec{8, 6, 4, 2};
    P(4, 3) = vec{4, 2, 0, 1};

    const auto A = U.evaluate_point(2.5, 1., P);

    REQUIRE(A(0) == Approx(2.));
    REQUIRE(A(1) == Approx(98. / 27.));
    REQUIRE(A(2) == Approx(68. / 27.));
}

TEST_CASE("NURBSSurface Compute Derivative of A Point", "[IGA.NURBSSurface]") {
    const auto U = NURBSSurface(vec{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}, vec{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}, 2);

    field<vec> P(8, 8);
    P(2, 2) = vec{2, 1};
    P(2, 3) = vec{6, 1};
    P(2, 4) = vec{7, 2};
    P(3, 2) = vec{1, 1};
    P(3, 3) = vec{9, 1};
    P(3, 4) = vec{3, 1};
    P(4, 2) = vec{4, 1};
    P(4, 3) = vec{5, 1};
    P(4, 4) = vec{8, 1};

    auto PP = P;

    IGA::convert_to_weighted(P);

    const auto A = U.evaluate_point_derivative(2.5, 2.5, P, 5, 5);

    REQUIRE(A(0, 0)(0) == Approx(34. / 5.));
    REQUIRE(A(1, 0)(0) == Approx(-64. / 325.));
    REQUIRE(A(2, 0)(0) == Approx(-75392. / 21125.));
    REQUIRE(A(3, 0)(0) == Approx(-804864. / 1373125.));
    REQUIRE(A(4, 0)(0) == Approx(222345216. / 89253125.));
    REQUIRE(A(1, 1)(0) == Approx(-4288. / 21125.));
    REQUIRE(A(0, 2)(0) == Approx(-223872. / 21125.));

    const auto B = U.evaluate_shape_function_derivative(2.5, 2.5, P, 2, 2)(2, 2);

    vec C(2, fill::zeros);

    for(auto I = 0; I < 64; ++I) if(!PP(I).empty()) C += PP(I) * B(I);

    REQUIRE(A(2, 2)(0) == Approx(C(0)));
}

TEST_CASE("NURBSVolume Compute Location of A Point", "[IGA.NURBSVolume]") {
    const auto U = NURBSVolume3D(vec{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}, vec{0, 0, 0, 1, 2, 3, 3, 3}, vec{0, 0, 0, 1, 2, 3, 3, 3});

    field<vec> P(10, 10, 10);
    P.for_each([](vec& N) { N.randu(4); });

    auto PP = P;

    IGA::convert_to_weighted(P);

    const auto A = U.evaluate_point(2.5, 1., .5, P);

    const auto B = U.evaluate_shape_function(2.5, 1., .5, P);

    vec C = zeros(4);

    for(auto I = 0llu; I < B.n_rows; ++I) for(auto J = 0llu; J < B.n_cols; ++J) for(auto K = 0llu; K < B.n_slices; ++K) C += B(I, J, K) * P(I, J, K);

    REQUIRE(norm(A) == Approx(norm(C.head(3))));

    const auto D = U.evaluate_point_derivative(2.5, 1., .5, P, 2)(2, 2, 2);

    const auto E = U.evaluate_shape_function_derivative(2.5, 1., .5, P, 2, 2, 2)(2, 2, 2);

    C.zeros();

    for(auto I = 0llu; I < E.n_rows; ++I) for(auto J = 0llu; J < E.n_cols; ++J) for(auto K = 0llu; K < E.n_slices; ++K) C += E(I, J, K) * PP(I, J, K);

    REQUIRE(norm(D) == Approx(norm(C.head(3))));
}
