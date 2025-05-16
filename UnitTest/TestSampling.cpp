#include "CatchHeader.h"

#include <Toolbox/resampling.h>

TEST_CASE("GCD", "[Utility.Sampling]") {
    REQUIRE(gcd(20, 10) == 10);
    REQUIRE(gcd(14, 24) == 2);
    REQUIRE(gcd(56, 98) == 14);
    REQUIRE(gcd(51248431, 3411571202) == 13);
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
