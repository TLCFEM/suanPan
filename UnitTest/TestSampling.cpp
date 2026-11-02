#include "CatchHeader.h"

#include <Toolbox/resampling.h>

TEST_CASE("GCD", "[Utility.Sampling]") {
    REQUIRE(gcd(20, 10) == 10);
    REQUIRE(gcd(14, 24) == 2);
    REQUIRE(gcd(56, 98) == 14);
    REQUIRE(gcd(51248431, 3411571202) == 13);
    REQUIRE(gcd(1, 1) == 1);
    REQUIRE(gcd(0, 5) == 5);
    REQUIRE(gcd(7, 0) == 7);
    REQUIRE(gcd(17, 13) == 1); // coprime
}

TEST_CASE("Window Functions", "[Utility.Sampling]") {
    // all windows evaluated at n=9 (indices 0..8) must produce symmetric, bounded values
    for(const auto& w : {hamming(9), hann(9), blackman(9), blackman_nuttall(9), blackman_harris(9), flat_top(9)}) {
        REQUIRE(9 == w.n_elem);
        // symmetry about the centre
        REQUIRE(Approx(w(0)).margin(1E-12) == w(8));
        REQUIRE(Approx(w(1)).margin(1E-12) == w(7));
        // bounds scaled to nominal range (with a small rounding allowance for the flat-top window)
        REQUIRE(w.min() >= -1. - 1E-6);
        REQUIRE(w.max() <= 1. + 1E-6);
    }

    // hamming endpoints: 0.54 - 0.46 = 0.08; centre (index 4 of 9): 0.54 + 0.46 = 1.0
    const auto hm = hamming(9);
    REQUIRE(Approx(.08).margin(1E-12) == hm(0));
    REQUIRE(Approx(1.).margin(1E-12) == hm(4));

    // hann endpoints: 0; centre: 1
    const auto hn = hann(9);
    REQUIRE(Approx(0.).margin(1E-12) == hn(0));
    REQUIRE(Approx(1.).margin(1E-12) == hn(4));
}

TEST_CASE("FIR Filter Normalised", "[Utility.Sampling]") {
    // low pass / band stop filters are normalised so the coefficients sum to 1 (DC gain 1)
    const auto lp = fir_low_pass<WindowType::Hamming>(8, .4);
    REQUIRE(Approx(1.).margin(1E-9) == sum(lp));

    const auto bs = fir_band_stop<WindowType::Hamming>(10, .2, .6);
    REQUIRE(Approx(1.).margin(1E-9) == sum(bs));

    // invalid order throws for even-order-only filters (only when assertions are active)
#ifdef SUANPAN_DEBUG
    REQUIRE_THROWS(fir_high_pass<WindowType::Hamming>(7, .4));
    REQUIRE_THROWS(fir_band_stop<WindowType::Hamming>(9, .2, .6));
    // invalid frequency ordering
    REQUIRE_THROWS(fir_band_pass<WindowType::Hamming>(8, .8, .1));
#endif
}

TEST_CASE("Upsampling Preserves Signal", "[Utility.Sampling]") {
    // a constant signal upsampled must remain constant (edge transients excluded)
    const vec in(32, fill::ones);
    vec out = upsampling<WindowType::Hamming>(in, 4, 8);
    REQUIRE(out.n_elem == 4 * in.n_elem);

    // examine the central 40% of the output, away from "same" convolution edge effects
    const auto lo = out.n_elem / 10 * 3;
    const auto hi = out.n_elem / 10 * 7;
    const auto mean_val = mean(out.subvec(lo, hi));
    REQUIRE(Approx(1.).margin(.05) == mean_val);

    const auto max_dev = max(abs(out.subvec(lo, hi)) - 1.);
    REQUIRE(max_dev < .15);
}

TEST_CASE("FIR", "[Utility.Sampling]") {
    fir_low_pass<WindowType::Hamming>(4, .1);
    fir_low_pass<WindowType::Hann>(4, .1);
    fir_low_pass<WindowType::Blackman>(4, .1);
    fir_low_pass<WindowType::BlackmanNuttall>(4, .1);
    fir_low_pass<WindowType::BlackmanHarris>(4, .1);
    fir_low_pass<WindowType::FlatTop>(4, .1);

    fir_high_pass<WindowType::Hamming>(4, .1);
    fir_high_pass<WindowType::Hann>(4, .1);
    fir_high_pass<WindowType::Blackman>(4, .1);
    fir_high_pass<WindowType::BlackmanNuttall>(4, .1);
    fir_high_pass<WindowType::BlackmanHarris>(4, .1);
    fir_high_pass<WindowType::FlatTop>(4, .1);

    fir_band_pass<WindowType::Hamming>(4, .1, .8);
    fir_band_pass<WindowType::Hann>(4, .1, .8);
    fir_band_pass<WindowType::Blackman>(4, .1, .8);
    fir_band_pass<WindowType::BlackmanNuttall>(4, .1, .8);
    fir_band_pass<WindowType::BlackmanHarris>(4, .1, .8);
    fir_band_pass<WindowType::FlatTop>(4, .1, .8);

    fir_band_stop<WindowType::Hamming>(4, .1, .8);
    fir_band_stop<WindowType::Hann>(4, .1, .8);
    fir_band_stop<WindowType::Blackman>(4, .1, .8);
    fir_band_stop<WindowType::BlackmanNuttall>(4, .1, .8);
    fir_band_stop<WindowType::BlackmanHarris>(4, .1, .8);
    fir_band_stop<WindowType::FlatTop>(4, .1, .8);
}
