/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
#ifdef __cpp_lib_execution
#include <execution>
#endif

namespace suanpan {
    template<sp_i IT, std::invocable<IT> F> void for_each(const IT start, const IT end, F&& FN) {
#ifdef SUANPAN_MT
        static tbb::affinity_partitioner ap;
        tbb::parallel_for(start, end, std::forward<F>(FN), ap);
#else
        for(IT I = start; I < end; ++I) FN(I);
#endif
    }

    template<sp_i IT, std::invocable<IT> F> void for_each(const IT end, F&& FN) { return for_each(static_cast<IT>(0), end, std::forward<F>(FN)); }

    template<typename T> constexpr T max_element(T start, T end) {
#ifdef __cpp_lib_execution
        return std::max_element(std::execution::par, start, end);
#else
        return std::max_element(start, end);
#endif
    }

    template<typename T> [[maybe_unused]] const std::vector<T>& unique(std::vector<T>& container) {
        std::sort(container.begin(), container.end());
        container.erase(std::unique(container.begin(), container.end()), container.end());
        container.shrink_to_fit();
        return container;
    }

    template<typename T> constexpr T& hacker(const T& I) { return const_cast<T&>(I); }

    template<typename T> constexpr T*& hacker(const T* const& I) { return const_cast<T*&>(I); }

    template<typename T> T sign(const T& I) { return (I > T(0)) - (I < T(0)); }

    template<typename T> bool approx_equal(T x, T y, int ulp = 2) requires (!std::numeric_limits<T>::is_integer) { return fabs(x - y) <= std::numeric_limits<T>::epsilon() * fabs(x + y) * ulp || fabs(x - y) < std::numeric_limits<T>::min(); }

    unsigned long long binomial(unsigned long long, unsigned long long);

    char to_upper(char);
    char to_lower(char);

    void to_upper(string&);
    void to_lower(string&);
    string to_upper(const string&);
    string to_lower(const string&);
    string to_upper(string&&);
    string to_lower(string&&);

    namespace expression {
        std::vector<std::pair<string, unsigned>> split(const std::string_view& variable_string);
    } // namespace expression 
}     // namespace suanpan

template<typename T> bool get_input(istringstream& I, T& O) { return static_cast<bool>(I >> O); }

template<typename T> bool get_input(istringstream& I, Col<T>& O) {
    auto code = true;
    for(auto& P : O) code &= static_cast<bool>(I >> P);
    return code;
}

template<typename T> bool get_input(istringstream& I, std::vector<T>& O) {
    T value;
    while(get_input(I, value)) O.emplace_back(value);
    return true;
}

template<typename T, typename... U> bool get_input(istringstream& I, T& O, U&... R) { return static_cast<bool>(I >> O) ? get_input(I, R...) : false; }

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

template<typename T, typename... U> bool get_optional_input(istringstream& I, T& O, U&... R) {
    if(I.eof()) return true;

    return static_cast<bool>(I >> O) ? get_optional_input(I, R...) : false;
}

string get_remaining(istringstream&);

bool is_equal(const char*, const char*);
bool is_equal(char, char);
bool is_equal(int, char);
bool is_equal(const string&, const char*);
bool is_equal(const char*, const string&);
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

bool is_integer(const string&);

double perturb(double, double = 1E-5);

#endif
