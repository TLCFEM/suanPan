#include "CatchHeader.h"

#include <Toolbox/utility.h>
#include <sstream>

using namespace suanpan;

TEST_CASE("Binomial Compute Basic Function", "[Utility.Binomial]") {
    REQUIRE(165 == suanpan::binomial(11, 3));
    REQUIRE(6435 == suanpan::binomial(15, 7));
    REQUIRE(455 == suanpan::binomial(15, 3));
    REQUIRE(126 == suanpan::binomial(9, 4));
    REQUIRE(220 == suanpan::binomial(12, 3));
}

TEST_CASE("Binomial Edge Cases", "[Utility.Binomial]") {
    REQUIRE(1 == suanpan::binomial(0, 0));
    REQUIRE(1 == suanpan::binomial(10, 0));
    REQUIRE(1 == suanpan::binomial(10, 10));
    REQUIRE(10 == suanpan::binomial(10, 1));
    // symmetry: C(n,k) == C(n,n-k)
    REQUIRE(suanpan::binomial(10, 2) == suanpan::binomial(10, 8));
    // pascal triangle row 5: 1 5 10 10 5 1
    REQUIRE(5 == suanpan::binomial(5, 1));
    REQUIRE(10 == suanpan::binomial(5, 2));
    REQUIRE(10 == suanpan::binomial(5, 3));
    REQUIRE(5 == suanpan::binomial(5, 4));
}

TEST_CASE("Sign", "[Utility.Sign]") {
    REQUIRE(Approx(1) == suanpan::sign(std::numeric_limits<double>::epsilon()));
    REQUIRE(Approx(-1) == suanpan::sign(-std::numeric_limits<double>::epsilon()));
    REQUIRE(Approx(0) == suanpan::sign(0));
    REQUIRE(1 == suanpan::sign(42));
    REQUIRE(-1 == suanpan::sign(-3));
    REQUIRE(1 == suanpan::sign(3.14));
    REQUIRE(-1 == suanpan::sign(-2.71));
    REQUIRE(Approx(0) == suanpan::sign(0.));
}

TEST_CASE("Clamp", "[Utility.Clamp]") {
    REQUIRE(Approx(2.) == suanpan::clamp(5., 0., 2.));
    REQUIRE(Approx(0.) == suanpan::clamp(-5., 0., 2.));
    REQUIRE(Approx(1.5) == suanpan::clamp(1.5, 0., 2.));
    // arguments may be out of order
    REQUIRE(Approx(2.) == suanpan::clamp(5., 2., 0.));
    REQUIRE(Approx(0.) == suanpan::clamp(-5., 2., 0.));
    // clamp with explicit test/bound ranges
    // value within (low_test, high_test) is passed through unchanged
    REQUIRE(Approx(.5) == suanpan::clamp(.5, .2, .8, 1., 2.));
    // value above high_test maps to high_bound
    REQUIRE(Approx(2.) == suanpan::clamp(.9, .2, .8, 1., 2.));
    // value below low_test maps to low_bound
    REQUIRE(Approx(1.) == suanpan::clamp(.1, .2, .8, 1., 2.));
}

TEST_CASE("Clamp Unit", "[Utility.Clamp]") {
    REQUIRE(Approx(0.) == suanpan::clamp_unit<double>(-.5));
    REQUIRE(Approx(1.) == suanpan::clamp_unit<double>(1.5));
    REQUIRE(Approx(.5) == suanpan::clamp_unit<double>(.5));
    REQUIRE(Approx(0.) == suanpan::clamp_unit(0.));
}

TEST_CASE("Approx Equal", "[Utility.Approx]") {
    REQUIRE(suanpan::approx_equal(1., 1. + 1E-16));
    REQUIRE(!suanpan::approx_equal(1., 1. + 1E-3));
    REQUIRE(suanpan::approx_equal(0., 0.));
    // within ulp tolerance for large magnitudes
    REQUIRE(suanpan::approx_equal(1000., 1000. + 1E-13, 100));
    REQUIRE(!suanpan::approx_equal(1., 2.));
}

TEST_CASE("Middle Element", "[Utility.Middle]") {
    REQUIRE(2 == suanpan::middle(std::vector<int>{1, 2, 3}));
    REQUIRE(3 == suanpan::middle(std::vector<int>{1, 2, 3, 4})); // element at size/2 index
    REQUIRE(5 == suanpan::middle(std::vector<double>{1., 5., 9.}));
}

TEST_CASE("String Case Conversion", "[Utility.StringCase]") {
    REQUIRE('A' == suanpan::to_upper('a'));
    REQUIRE('A' == suanpan::to_upper('A'));
    REQUIRE('z' == suanpan::to_lower('Z'));
    REQUIRE('z' == suanpan::to_lower('z'));

    std::string s = "AbC dEf";
    to_upper(s);
    REQUIRE("ABC DEF" == s);
    to_lower(s);
    REQUIRE("abc def" == s);

    REQUIRE("ABC" == suanpan::to_upper(std::string("abc")));
    REQUIRE("abc" == suanpan::to_lower(std::string("ABC")));
}

TEST_CASE("Unique Container", "[Utility.Unique]") {
    auto v = std::vector{3, 1, 2, 2, 3, 1, 4};
    const auto& r = suanpan::unique(v);
    REQUIRE(v.size() == 4);
    REQUIRE(r.size() == 4);
    REQUIRE(v[0] == 1);
    REQUIRE(v[1] == 2);
    REQUIRE(v[2] == 3);
    REQUIRE(v[3] == 4);
}

TEST_CASE("Max Element", "[Utility.MaxElement]") {
    std::vector v{1, 5, 3, 9, 2};
    REQUIRE(9 == *suanpan::max_element(v.begin(), v.end()));
    std::vector w{1.1, .5, 2.2};
    REQUIRE(Approx(2.2) == *suanpan::max_element(w.begin(), w.end()));
}

TEST_CASE("Parallel All/None/Any Of", "[Utility.Parallel]") {
    std::vector v{2, 4, 6, 8};
    REQUIRE(suanpan::all_of(v, [](const int i) { return i % 2 == 0; }));
    REQUIRE(!suanpan::all_of(v, [](const int i) { return i > 4; }));
    REQUIRE(suanpan::none_of(v, [](const int i) { return i % 2 == 1; }));
    REQUIRE(!suanpan::none_of(v, [](const int i) { return i > 4; }));
    REQUIRE(suanpan::any_of(v, [](const int i) { return i == 8; }));
    REQUIRE(!suanpan::any_of(v, [](const int i) { return i == 9; }));
}

TEST_CASE("For Each", "[Utility.ForEach]") {
    std::vector v(10, 0);
    for_each(static_cast<uword>(v.size()), [&](const uword I) { ++v[I]; });
    for(const auto X : v) REQUIRE(1 == X);

    std::vector w(5, 0);
    for_each(1u, 5u, [&](const uword I) { ++w[I]; });
    REQUIRE(0 == w[0]);
    for(auto I = 1u; I < 5u; ++I) REQUIRE(1 == w[I]);
}

TEST_CASE("Hacker Const Cast", "[Utility.Hacker]") {
    const auto value = 42;
    const auto& cref = value;
    hacker(cref) = 43;
    REQUIRE(43 == value);

    const auto ptr = new int(3);
    const auto& cpref = ptr;
    const auto original = ptr;
    hacker(cpref) = new int(5);
    REQUIRE(5 == *ptr);
    delete original;
    delete ptr;
}

TEST_CASE("Get Input Scalar", "[Utility.GetInput]") {
    std::istringstream stream("1 2 3");
    double a{}, b{}, c{};
    REQUIRE(get_input(stream, a, b, c));
    REQUIRE(Approx(1) == a);
    REQUIRE(Approx(2) == b);
    REQUIRE(Approx(3) == c);
}

TEST_CASE("Get Input Vector", "[Utility.GetInput]") {
    std::istringstream stream("1 2 3");
    vec out(3, fill::none);
    REQUIRE(get_input(stream, out));
    REQUIRE(Approx(1) == out(0));
    REQUIRE(Approx(2) == out(1));
    REQUIRE(Approx(3) == out(2));
}

TEST_CASE("Get Input Tuple", "[Utility.GetInput]") {
    std::istringstream stream("1 2.5 x");
    std::tuple<int, double, std::string> out;
    REQUIRE(get_input(stream, out));
    REQUIRE(1 == std::get<0>(out));
    REQUIRE(Approx(2.5) == std::get<1>(out));
    REQUIRE("x" == std::get<2>(out));
}

TEST_CASE("Get Remaining Vector", "[Utility.GetRemaining]") {
    std::istringstream stream("1 2 3 4");
    const auto v = get_remaining<int>(stream);
    REQUIRE(4 == v.size());
    REQUIRE(1 == v[0]);
    REQUIRE(4 == v[3]);
}

TEST_CASE("Get Remaining As Set", "[Utility.GetRemainingSet]") {
    std::istringstream stream("3 1 2 3 1");
    const auto s = get_remaining_as_set<int>(stream);
    REQUIRE(3 == s.size());
    REQUIRE(1 == *s.begin());
    REQUIRE(3 == *s.rbegin());
}

TEST_CASE("get_input Templated Stream", "[Utility.GetInputT]") {
    std::istringstream stream("7");
    REQUIRE(7 == get_input<int>(stream));
}

TEST_CASE("Get Remaining Multi Type", "[Utility.GetRemainingMulti]") {
    // number of columns in each vector must match
    std::istringstream stream("1 2\n3 4\n5 6");
    const auto [ints, doubles] = get_remaining<int, double>(stream);
    REQUIRE(3 == ints.size());
    REQUIRE(3 == doubles.size());
    REQUIRE(1 == ints[0]);
    REQUIRE(4 == doubles[1]);
}

TEST_CASE("Get Remaining As Tuple", "[Utility.GetRemainingTuple]") {
    std::istringstream stream("1 2\n3 4");
    const auto v = get_remaining_as_tuple<int, double>(stream);
    REQUIRE(2 == v.size());
    REQUIRE(3 == std::get<0>(v[1]));
    REQUIRE(Approx(4) == std::get<1>(v[1]));
}

TEST_CASE("Ignore Whitespace", "[Utility.IgnoreWhitespace]") {
    std::istringstream stream("   \t\tabc");
    ignore_whitespace(stream);
    REQUIRE('a' == stream.peek());
}

TEST_CASE("Get Remaining String", "[Utility.GetRemainingString]") {
    std::istringstream stream("  hello world ");
    const auto s = get_remaining(stream);
    REQUIRE("hello world " == s);
}

TEST_CASE("Get Optional Input", "[Utility.GetOptional]") {
    std::istringstream stream("1 2");
    double a{}, b{}, c{};
    REQUIRE(get_optional_input(stream, a, b, c));
    REQUIRE(Approx(1) == a);
    REQUIRE(Approx(2) == b);
    REQUIRE(Approx(0) == c);

    std::istringstream stream2("5");
    double d{};
    REQUIRE(get_optional_input(stream2, d));
    REQUIRE(Approx(5) == d);
    REQUIRE(stream2.eof());
}

TEST_CASE("Is Equal", "[Utility.IsEqual]") {
    REQUIRE(is_equal('a', 'A'));
    REQUIRE(is_equal('a', 'a'));
    REQUIRE(!is_equal('a', 'b'));
    // int-vs-char overload: numeric code of 'a' (97) equals 'A' case-insensitively
    REQUIRE(is_equal(static_cast<int>('a'), 'A'));
    REQUIRE(!is_equal(static_cast<int>('a'), 'B'));
    REQUIRE(is_equal(std::string_view("Hello"), std::string_view("hELLo")));
    REQUIRE(!is_equal(std::string_view("Hello"), std::string_view("World")));
}

TEST_CASE("Is Equal Any", "[Utility.IsEqualAny]") {
    REQUIRE(is_equal_any(std::string_view("X"), "x", "y"));
    REQUIRE(!is_equal_any(std::string_view("Z"), "x", "y"));
}

TEST_CASE("If Startswith", "[Utility.IfStartswith]") {
    REQUIRE(if_startswith("suanPan", "suan"));
    REQUIRE(if_startswith("suanPan", "SUAN"));
    REQUIRE(!if_startswith("suanPan", "pan"));
    REQUIRE(!if_startswith("ab", "abc"));
}

TEST_CASE("If Contain", "[Utility.IfContain]") {
    REQUIRE(if_contain(std::string("hello world"), "world"));
    REQUIRE(if_contain(std::string("hello world"), std::string("lo wo")));
    REQUIRE(!if_contain(std::string("hello"), "world"));
}

TEST_CASE("Is True / Is False", "[Utility.IsTrue]") {
    for(const std::vector<std::string> truthy{"On", "True", "T", "1", "Yes", "Y", "yes", "tRuE"}; const auto& s : truthy) REQUIRE(is_true(s));

    for(const std::vector<std::string> falsy{"Off", "False", "F", "0", "No", "N", "no", "fAlSe"}; const auto& s : falsy) REQUIRE(is_false(s));

    REQUIRE(!is_true("Maybe"));
    REQUIRE(!is_false("Maybe"));
}

TEST_CASE("Is Integer", "[Utility.IsInteger]") {
    REQUIRE(is_integer("0"));
    REQUIRE(is_integer("12345"));
    REQUIRE(!is_integer(""));
    REQUIRE(!is_integer("12a"));
    REQUIRE(!is_integer("1.5"));
    REQUIRE(!is_integer("-1"));
}

TEST_CASE("Perturb", "[Utility.Perturb]") {
    // perturb scales by 1 + randu * fraction, so the result must remain within the expected band
    for(auto I = 0; I < 100; ++I) {
        const auto p = perturb(10., .01);
        REQUIRE(p >= 10.);
        REQUIRE(p <= 10. + 10. * .01);
        REQUIRE(Approx(10.).margin(0.11) == p);
    }
}

TEST_CASE("Expression Split Empty / Invalid", "[Utility.ExpressionSplit]") {
    const auto r1 = expression::split("x|y");
    REQUIRE(2 == r1.size());
    REQUIRE(1u == r1[0].second);
    REQUIRE(1u == r1[1].second);

    // leading integer means invalid expression
    const auto r2 = expression::split("5|x");
    REQUIRE(r2.empty());
}

TEST_CASE("Matrix Allocation", "[Utility.Matrix]") {
    BENCHMARK("Static Size 20") {
        mat::fixed<20, 20> A(fill::randn);
        A(10, 10) = 1.;
        REQUIRE(A.n_elem == 400);
        return A;
    };

    BENCHMARK("Dynamic Size 20") {
        mat A(20, 20, fill::randn);
        A(10, 10) = 1.;
        REQUIRE(A.n_elem == 400);
        return A;
    };
}

TEST_CASE("Matrix Resize", "[Utility.Matrix]") {
    BENCHMARK("Upsize") { return vec(4).resize(6); };

    BENCHMARK("Downsize") { return vec(4).resize(2); };
}

TEST_CASE("Color Print", "[Utility.Print]") {
    suanpan_info("TEST.\n");
    suanpan_highlight("TEST.\n");
    suanpan_debug("TEST.\n");
    suanpan_warning("TEST.\n");
    suanpan_error("TEST.\n");
    suanpan_fatal("TEST.\n");
    suanpan_info("TEST.\n", vec{1, 2, 3});
}
