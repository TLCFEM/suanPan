/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#ifndef RESAMPLING_H
#define RESAMPLING_H

#include <Toolbox/utility.h>

enum class WindowType {
    Hamming,
    Hann,
    Blackman,
    BlackmanNuttall,
    BlackmanHarris,
    FlatTop
};

uword gcd(uword, uword);

vec hamming(uword);
vec hann(uword);
vec blackman(uword);
vec blackman_nuttall(uword);
vec blackman_harris(uword);
vec flat_top(uword);

vec fir_low_pass(uword, double, vec (*)(uword));

template<WindowType T> vec fir_low_pass(uword, double) { throw invalid_argument("unknown window type"); }

template<> vec fir_low_pass<WindowType::Hamming>(uword, double);

template<> vec fir_low_pass<WindowType::Hann>(uword, double);

template<> vec fir_low_pass<WindowType::Blackman>(uword, double);

template<> vec fir_low_pass<WindowType::BlackmanNuttall>(uword, double);

template<> vec fir_low_pass<WindowType::BlackmanHarris>(uword, double);

template<> vec fir_low_pass<WindowType::FlatTop>(uword, double);

template<WindowType T> vec upsampling(const vec& in, const uword up_rate) {
    vec out(up_rate * in.n_elem, fill::none);

    const auto coef = fir_low_pass<T>(8llu * up_rate, 1. / static_cast<double>(up_rate));

    const auto buffer_size = coef.n_elem;

    vec buffer(buffer_size, fill::zeros);
    auto index = 0llu;

    auto feed = [&](const double next) {
        buffer(index) = next;

        uword p = 0;
        double result = 0.;
        for(auto m = index; m < buffer_size; ++m) result += coef(p++) * buffer(m);
        for(auto m = 0llu; m < index; ++m) result += coef(p++) * buffer(m);

        index = (index == 0 ? buffer_size : index) - 1;

        return result;
    };

    for(auto n = 0llu; n < out.n_elem; ++n) out(n) = static_cast<double>(up_rate) * feed(0 == n % up_rate ? in(n / up_rate) : 0.);

    return out;
}

template<WindowType T> mat upsampling(const string& file_name, const uword up_rate) {
    mat ext_data;

    if(!fs::exists(file_name) || !ext_data.load(file_name, raw_ascii) || ext_data.empty() || ext_data.n_cols < 2) {
        ext_data.reset();
        return ext_data;
    }

    const vec time_diff = diff(ext_data.col(0));

    if(!suanpan::approx_equal(time_diff.min(), time_diff.max(), 1000000)) {
        ext_data.reset();
        return ext_data;
    }

    mat result(ext_data.n_rows * up_rate, 2, fill::none);

    result.col(1) = upsampling<T>(ext_data.col(1), up_rate);

    const auto time_size = mean(time_diff) / static_cast<double>(up_rate);

    for(auto I = 0llu; I < result.n_rows; ++I) result(I, 0) = static_cast<double>(I) * time_size;

    return result;
}

#endif
