#include "CatchHeader.h"
#include <Section/SectionNM/SurfaceNM2D.h>

TEST_CASE("2D Surface", "[NM.Surface]") {
    const SurfaceNM2D A(1.);

    REQUIRE(Approx(.8406756) == A.compute_sf(vec{.5, .9}, 1.));

    const auto DF = A.compute_dsf(vec{.5, .9}, 1.);

    REQUIRE(Approx(4.122700) == DF(0));
    REQUIRE(Approx(3.451500) == DF(1));

    const auto DDF = A.compute_ddsf(vec{.5, .9}, 1.);

    REQUIRE(Approx(8.245400) == DDF(0, 0));
    REQUIRE(Approx(3.835000) == DDF(1, 1));
    REQUIRE(Approx(6.606000) == DDF(0, 1));
    REQUIRE(Approx(6.606000) == DDF(1, 0));
}
