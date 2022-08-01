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

#ifndef UTILITY_H
#define UTILITY_H

#include <suanPan.h>
#include <concepts>

template<sp_i IT, typename F> void suanpan_for(const IT start, const IT end, F&& FN) {
#ifdef SUANPAN_MT
    tbb::parallel_for(start, end, std::forward<F>(FN));
#else
    for(IT I = start; I < end; ++I) FN(I);
#endif
}

namespace suanpan {
    template<typename T> [[maybe_unused]] const std::vector<T>& unique(std::vector<T>& container) {
        std::sort(container.begin(), container.end());
        container.erase(std::unique(container.begin(), container.end()), container.end());
        container.shrink_to_fit();
        return container;
    }

    template<typename T> constexpr T& hacker(const T& I) { return const_cast<T&>(I); }

    template<typename T> constexpr T*& hacker(const T* const& I) { return const_cast<T*&>(I); }

    template<typename T> T sign(const T& I) { return (I > T(0)) - (I < T(0)); }

    template<typename T> std::enable_if_t<!std::numeric_limits<T>::is_integer, bool> approx_equal(T x, T y, int ulp = 2) { return fabs(x - y) <= std::numeric_limits<T>::epsilon() * fabs(x + y) * ulp || fabs(x - y) < std::numeric_limits<T>::min(); }

    unsigned long long binomial(unsigned long long, unsigned long long);

    char to_upper(char);
    char to_lower(char);

    void to_upper(string&);
    void to_lower(string&);
    string to_upper(const string&);
    string to_lower(const string&);
    string to_upper(string&&);
    string to_lower(string&&);
} // namespace suanpan

template<typename T> bool get_input(istringstream& I, T& O) { return static_cast<bool>(I >> O); }

template<typename T> bool get_input(istringstream& I, Col<T>& O) {
    auto code = true;
    for(auto& P : O) code &= static_cast<bool>(I >> P);
    return code;
}

template<typename T, typename...U> bool get_input(istringstream& I, T& O, U&...R) { return static_cast<bool>(I >> O) ? get_input(I, R...) : false; }

template<typename T> T get_input(istringstream& I) {
    T O;
    I >> O;
    return O;
}

void ignore_whitespace(istringstream&);

template<typename T> bool get_optional_input(istringstream& I, T& O) {
    if(I.eof()) return true;

    return static_cast<bool>(I >> O);
}

template<typename T> bool get_optional_input(istringstream& I, Col<T>& O) {
    auto code = true;
    for(auto& P : O) code &= I.eof() ? true : static_cast<bool>(I >> P);
    return code;
}

template<typename T, typename...U> bool get_optional_input(istringstream& I, T& O, U&...R) {
    if(I.eof()) return true;

    return static_cast<bool>(I >> O) ? get_optional_input(I, R...) : false;
}

bool is_equal(const char*, const char*);
bool is_equal(char, char);
bool is_equal(int, char);
bool is_equal(const string&, const char*);
bool is_equal(const string&, const string&);

bool if_contain(const string&, const char*);
bool if_contain(const string&, const string&);
bool if_contain(string&&, string&&);

template<std::equality_comparable T> std::pair<bool, std::int64_t> if_contain(const std::vector<T>& container, const T target) {
    auto position = std::find(container.begin(), container.end(), target);

    return {position != container.end() && container.size() > 0, position - container.begin()};
}

bool is_true(const char*);
bool is_false(const char*);
bool is_true(const string&);
bool is_false(const string&);

#endif
