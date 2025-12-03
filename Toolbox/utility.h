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

#include <concepts>
#include <suanPan.h>
#ifdef __cpp_lib_execution
#include <execution>
#endif

namespace suanpan {
    template<typename T> constexpr const T& middle(const std::vector<T>& container) { return container[container.size() / 2]; }

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

    template<class Container, class Handler> requires requires(Container x) { x.begin(); x.end(); } auto all_of(Container& target, Handler&& func) {
        return std::all_of(
#ifdef __cpp_lib_execution
            std::execution::par,
#endif
            target.begin(), target.end(), std::forward<Handler>(func)
        );
    }

    template<class Container, class Handler> requires requires(Container x) { x.begin(); x.end(); } auto none_of(Container& target, Handler&& func) {
        return std::none_of(
#ifdef __cpp_lib_execution
            std::execution::par,
#endif
            target.begin(), target.end(), std::forward<Handler>(func)
        );
    }

    template<class Container, class Handler> requires requires(Container x) { x.begin(); x.end(); } auto any_of(Container& target, Handler&& func) {
        return std::any_of(
#ifdef __cpp_lib_execution
            std::execution::par,
#endif
            target.begin(), target.end(), std::forward<Handler>(func)
        );
    }

    template<typename T> constexpr T& hacker(const T& I) { return const_cast<T&>(I); }

    template<typename T> constexpr T*& hacker(const T* const& I) { return const_cast<T*&>(I); }

    template<typename T> requires std::signed_integral<T> || std::floating_point<T> constexpr T sign(const T I) { return (I > T(0)) - (I < T(0)); }

    template<std::floating_point T> constexpr T clamp(const T c, T a, T b) {
        if(a > b) std::swap(a, b);
        return std::max(a, std::min(b, c));
    }

    template<std::floating_point T> constexpr T clamp_unit(const T c) { return clamp(c, T(0), T(1)); }

    template<std::floating_point T> bool approx_equal(const T x, const T y, int ulp = 2) { return std::fabs(x - y) <= std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp || std::fabs(x - y) < std::numeric_limits<T>::min(); }

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
} // namespace suanpan

template<typename T> bool get_input(std::istringstream& stream, T& output) { return static_cast<bool>(stream >> output); }

template<typename T> bool get_input(std::istringstream& stream, Col<T>& output) {
    for(auto& item : output)
        if(!get_input(stream, item)) return false;
    return true;
}

template<typename T, typename... U> bool get_input(std::istringstream& stream, T& output, U&... rest) { return get_input(stream, output) && get_input(stream, rest...); }

template<typename... Ts> bool get_input(std::istringstream& stream, std::tuple<Ts...>& output) {
    return std::apply([&](auto&... items) { return get_input(stream, items...); }, output);
}

template<typename T> auto get_remaining(std::istringstream& stream) {
    std::vector<T> output;
    T value;
    while(get_input(stream, value)) output.emplace_back(value);
    return output;
}

template<typename T1, typename T2, typename... Ts> auto get_remaining(std::istringstream& stream) {
    std::tuple<std::vector<T1>, std::vector<T2>, std::vector<Ts>...> output;
    std::tuple<T1, T2, Ts...> temp;
    while(get_input(stream, temp)) std::apply([&](auto&... containers) { std::apply([&](auto&&... items) { ((containers.emplace_back(items)), ...); }, temp); }, output);
    return output;
}

template<typename... Ts> auto get_remaining_as_tuple(std::istringstream& stream) {
    std::vector<std::tuple<Ts...>> output;
    std::tuple<Ts...> temp;
    while(get_input(stream, temp)) output.emplace_back(temp);
    return output;
}

template<typename T> T get_input(std::istringstream& stream) {
    T output{};
    stream >> output;
    return output;
}

void ignore_whitespace(std::istringstream&);

template<typename T> bool get_optional_input(std::istringstream& stream, T& output) { return stream.eof() || get_input(stream, output); }

template<typename T> bool get_optional_input(std::istringstream& stream, Col<T>& output) {
    for(auto& item : output)
        if(!stream.eof() && !get_input(stream, item)) return false;
    return true;
}

template<typename T, typename... U> bool get_optional_input(std::istringstream& stream, T& output, U&... rest) { return stream.eof() || (get_input(stream, output) && get_optional_input(stream, rest...)); }

std::string get_remaining(std::istringstream&);

bool is_equal(char, char);
bool is_equal(int, char);
bool is_equal(std::string_view, std::string_view);

template<typename... S> bool is_equal_any(std::string_view a, S... rest) { return (is_equal(a, rest) || ...); }

bool if_contain(const std::string&, const char*);
bool if_contain(const std::string&, const std::string&);

bool is_true(std::string_view);
bool is_false(std::string_view);

bool is_integer(const std::string&);

double perturb(double, double = 1E-5);

#endif
