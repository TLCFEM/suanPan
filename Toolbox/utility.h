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

    void to_upper(std::string&);
    void to_lower(std::string&);
    std::string to_upper(const std::string&);
    std::string to_lower(const std::string&);
    std::string to_upper(std::string&&);
    std::string to_lower(std::string&&);

    namespace expression {
        std::vector<std::pair<std::string, unsigned>> split(std::string_view variable_string);
    } // namespace expression 
}     // namespace suanpan

template<typename T> bool get_input(istringstream& I, T& O) { return static_cast<bool>(I >> O); }

template<typename T> bool get_input(istringstream& I, Col<T>& O) {
    auto code = true;
    for(auto& P : O) code &= static_cast<bool>(I >> P);
    return code;
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

template<typename T> auto get_remaining(istringstream& I) {
    std::vector<T> O;
    T value;
    while(get_input(I, value)) O.emplace_back(value);
    return O;
}

template<typename T1, typename T2> auto get_remaining(istringstream& I) {
    std::vector<T1> O1;
    std::vector<T2> O2;
    T1 V1;
    T2 V2;
    while(get_input(I, V1, V2)) {
        O1.emplace_back(V1);
        O2.emplace_back(V2);
    }
    return std::make_tuple(O1, O2);
}

template<typename T1, typename T2, typename T3> auto get_remaining(istringstream& I) {
    std::vector<T1> O1;
    std::vector<T2> O2;
    std::vector<T3> O3;
    T1 V1;
    T2 V2;
    T3 V3;
    while(get_input(I, V1, V2, V3)) {
        O1.emplace_back(V1);
        O2.emplace_back(V2);
        O3.emplace_back(V3);
    }
    return std::make_tuple(O1, O2, O3);
}

template<typename T1, typename T2, typename T3, typename T4> auto get_remaining(istringstream& I) {
    std::vector<T1> O1;
    std::vector<T2> O2;
    std::vector<T3> O3;
    std::vector<T4> O4;
    T1 V1;
    T2 V2;
    T3 V3;
    T4 V4;
    while(get_input(I, V1, V2, V3, V4)) {
        O1.emplace_back(V1);
        O2.emplace_back(V2);
        O3.emplace_back(V3);
        O4.emplace_back(V4);
    }
    return std::make_tuple(O1, O2, O3, O4);
}

template<typename T1, typename T2, typename T3, typename T4, typename T5> auto get_remaining(istringstream& I) {
    std::vector<T1> O1;
    std::vector<T2> O2;
    std::vector<T3> O3;
    std::vector<T4> O4;
    std::vector<T5> O5;
    T1 V1;
    T2 V2;
    T3 V3;
    T4 V4;
    T5 V5;
    while(get_input(I, V1, V2, V3, V4, V5)) {
        O1.emplace_back(V1);
        O2.emplace_back(V2);
        O3.emplace_back(V3);
        O4.emplace_back(V4);
        O5.emplace_back(V5);
    }
    return std::make_tuple(O1, O2, O3, O4, O5);
}

std::string get_remaining(istringstream&);

bool is_equal(const char*, const char*);
bool is_equal(char, char);
bool is_equal(int, char);
bool is_equal(const std::string&, const char*);
bool is_equal(const char*, const std::string&);
bool is_equal(const std::string&, const std::string&);
bool is_equal(std::string_view, const char*);
bool is_equal(const char*, std::string_view);

bool if_contain(const std::string&, const char*);
bool if_contain(const std::string&, const std::string&);
bool if_contain(std::string&&, std::string&&);

template<std::equality_comparable T> std::pair<bool, std::int64_t> if_contain(const std::vector<T>& container, const T target) {
    auto position = std::find(container.begin(), container.end(), target);

    return {position != container.end() && container.size() > 0, position - container.begin()};
}

bool is_true(const char*);
bool is_false(const char*);
bool is_true(const std::string&);
bool is_false(const std::string&);

bool is_integer(const std::string&);

double perturb(double, double = 1E-5);

#endif
