/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "resampling.h"

uword gcd(uword a, uword b) {
    if(b > a) std::swap(a, b);

    while(b != 0llu) {
        const auto c = a % b;
        a = b;
        b = c;
    }

    return a;
}

vec cos_window(const uword n, const vec& a) {
    vec h(n);

    for(auto i = 0llu; i < n; ++i) {
        const auto term = datum::tau * static_cast<double>(i) / (static_cast<double>(n) - 1.);

        h(i) = a(0) - a(1) * std::cos(term) + a(2) * std::cos(2. * term) - a(3) * std::cos(3. * term) + a(4) * std::cos(4. * term);
    }

    return h;
}

vec hamming(const uword n) { return cos_window(n, {.54, .46, .0, .0, .0}); }

vec hann(const uword n) { return cos_window(n, {.5, .5, .0, .0, .0}); }

vec blackman(const uword n) { return cos_window(n, {7938. / 18608., 9240. / 18608., 1430. / 18608., .0, .0}); }

vec blackman_nuttall(const uword n) { return cos_window(n, {.3635819, .4891775, .1365995, .0106411, .0}); }

vec blackman_harris(const uword n) { return cos_window(n, {.35875, .48829, .14128, .01168, .0}); }

vec flat_top(const uword n) { return cos_window(n, {.21557895, .41663158, .277263158, .083578947, .006947368}); }

vec fir_low_pass(const uword s, const double f, vec (*window)(uword)) {
    const auto sp = s + 1;

    const auto h = window(sp);

    vec b(sp);

    for(auto m = 0llu; m < sp; ++m) b(m) = f * h(m) * sinc(f * (static_cast<double>(m) - .5 * static_cast<double>(s)));

    return b / sum(b);
}

vec fir_high_pass(const uword s, const double f, vec (*window)(uword)) {
    suanpan_assert([&] { if(0 != s % 2) throw invalid_argument("order must be even"); });

    const auto sp = s + 1;

    const auto h = window(sp);

    vec b(sp, fill::none);

    for(auto m = 0llu; m < sp; ++m) {
        const auto term = static_cast<double>(m) - .5 * static_cast<double>(s);
        b(m) = h(m) * (sinc(term) - f * sinc(f * term));
    }

    vec bb(sp);
    for(auto i = 0llu; i < sp; i += 2llu) bb(i) = b(i);
    for(auto i = 1llu; i < sp; i += 2llu) bb(i) = -b(i);

    return b / std::abs(sum(bb));
}

vec fir_band_pass(const uword s, const double fa, const double fb, vec (*window)(uword)) {
    suanpan_assert([&] { if(fb <= fa) throw invalid_argument("frequencies must be [0 < fa < fb < 1]"); });

    const auto sp = s + 1;

    const auto h = window(sp);

    vec b(sp, fill::none);

    for(auto m = 0llu; m < sp; ++m) {
        const auto term = static_cast<double>(m) - .5 * static_cast<double>(s);
        b(m) = h(m) * (fb * sinc(fb * term) - fa * sinc(fa * term));
    }

    const auto fc = .5 * (fa + fb);
    const std::complex pi(0., -datum::pi);
    const vec fv = regspace(0, static_cast<double>(s));

    return b / std::abs(sum(exp(pi * fv * fc) % b));
}

vec fir_band_stop(const uword s, const double fa, const double fb, vec (*window)(uword)) {
    suanpan_assert([&] {
        if(0 != s % 2) throw invalid_argument("order must be even");
        if(fb <= fa) throw invalid_argument("frequencies must be [0 < fa < fb < 1]");
    });

    const auto sp = s + 1;

    const auto h = window(sp);

    vec b(sp, fill::none);

    for(auto m = 0llu; m < sp; ++m) {
        const auto term = static_cast<double>(m) - .5 * static_cast<double>(s);
        b(m) = h(m) * (sinc(term) - fb * sinc(fb * term) + fa * sinc(fa * term));
    }

    return b / sum(b);
}

template<> vec fir_low_pass<WindowType::Hamming>(const uword s, const double f) { return fir_low_pass(s, f, hamming); }

template<> vec fir_low_pass<WindowType::Hann>(const uword s, const double f) { return fir_low_pass(s, f, hann); }

template<> vec fir_low_pass<WindowType::Blackman>(const uword s, const double f) { return fir_low_pass(s, f, blackman); }

template<> vec fir_low_pass<WindowType::BlackmanNuttall>(const uword s, const double f) { return fir_low_pass(s, f, blackman_nuttall); }

template<> vec fir_low_pass<WindowType::BlackmanHarris>(const uword s, const double f) { return fir_low_pass(s, f, blackman_harris); }

template<> vec fir_low_pass<WindowType::FlatTop>(const uword s, const double f) { return fir_low_pass(s, f, flat_top); }

template<> vec fir_high_pass<WindowType::Hamming>(const uword s, const double f) { return fir_high_pass(s, f, hamming); }

template<> vec fir_high_pass<WindowType::Hann>(const uword s, const double f) { return fir_high_pass(s, f, hann); }

template<> vec fir_high_pass<WindowType::Blackman>(const uword s, const double f) { return fir_high_pass(s, f, blackman); }

template<> vec fir_high_pass<WindowType::BlackmanNuttall>(const uword s, const double f) { return fir_high_pass(s, f, blackman_nuttall); }

template<> vec fir_high_pass<WindowType::BlackmanHarris>(const uword s, const double f) { return fir_high_pass(s, f, blackman_harris); }

template<> vec fir_high_pass<WindowType::FlatTop>(const uword s, const double f) { return fir_high_pass(s, f, flat_top); }

template<> vec fir_band_pass<WindowType::Hamming>(const uword s, const double fa, const double fb) { return fir_band_pass(s, fa, fb, hamming); }

template<> vec fir_band_pass<WindowType::Hann>(const uword s, const double fa, const double fb) { return fir_band_pass(s, fa, fb, hann); }

template<> vec fir_band_pass<WindowType::Blackman>(const uword s, const double fa, const double fb) { return fir_band_pass(s, fa, fb, blackman); }

template<> vec fir_band_pass<WindowType::BlackmanNuttall>(const uword s, const double fa, const double fb) { return fir_band_pass(s, fa, fb, blackman_nuttall); }

template<> vec fir_band_pass<WindowType::BlackmanHarris>(const uword s, const double fa, const double fb) { return fir_band_pass(s, fa, fb, blackman_harris); }

template<> vec fir_band_pass<WindowType::FlatTop>(const uword s, const double fa, const double fb) { return fir_band_pass(s, fa, fb, flat_top); }

template<> vec fir_band_stop<WindowType::Hamming>(const uword s, const double fa, const double fb) { return fir_band_stop(s, fa, fb, hamming); }

template<> vec fir_band_stop<WindowType::Hann>(const uword s, const double fa, const double fb) { return fir_band_stop(s, fa, fb, hann); }

template<> vec fir_band_stop<WindowType::Blackman>(const uword s, const double fa, const double fb) { return fir_band_stop(s, fa, fb, blackman); }

template<> vec fir_band_stop<WindowType::BlackmanNuttall>(const uword s, const double fa, const double fb) { return fir_band_stop(s, fa, fb, blackman_nuttall); }

template<> vec fir_band_stop<WindowType::BlackmanHarris>(const uword s, const double fa, const double fb) { return fir_band_stop(s, fa, fb, blackman_harris); }

template<> vec fir_band_stop<WindowType::FlatTop>(const uword s, const double fa, const double fb) { return fir_band_stop(s, fa, fb, flat_top); }

mat upsampling(const string& window_type, const string& file_name, const uword up_rate, const uword window_size) {
    mat result;

    if(is_equal(window_type, "Hamming")) result = upsampling<WindowType::Hamming>(file_name, up_rate, window_size);
    else if(is_equal(window_type, "Hann")) result = upsampling<WindowType::Hann>(file_name, up_rate, window_size);
    else if(is_equal(window_type, "Blackman")) result = upsampling<WindowType::Blackman>(file_name, up_rate, window_size);
    else if(is_equal(window_type, "BlackmanNuttall")) result = upsampling<WindowType::BlackmanNuttall>(file_name, up_rate, window_size);
    else if(is_equal(window_type, "BlackmanHarris")) result = upsampling<WindowType::BlackmanHarris>(file_name, up_rate, window_size);
    else if(is_equal(window_type, "FlatTop")) result = upsampling<WindowType::FlatTop>(file_name, up_rate, window_size);

    return result;
}
