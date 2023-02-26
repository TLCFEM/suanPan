/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#ifndef SUANPAN_H
#define SUANPAN_H

// SUANPAN_DEBUG
// _DEBUG --> MSVC
// DEBUG --> GCC
#if defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
#define SUANPAN_DEBUG
#define SUANPAN_EXTRA_DEBUG
#else
#define ARMA_NO_DEBUG
#endif

#ifdef SUANPAN_SUPERLUMT
#define ARMA_DONT_USE_SUPERLU
#else
#define ARMA_USE_SUPERLU
#endif

#ifdef SUANPAN_MKL
#define MKL_DIRECT_CALL
#endif

#ifdef SUANPAN_HDF5
#define ARMA_USE_HDF5
#endif

#ifdef SUANPAN_MAGMA
#ifndef SUANPAN_CUDA
#undef SUANPAN_MAGMA
#endif
#endif

// SUANPAN_WIN
// WIN32 _WIN32 __WIN32 __WIN32__ --> MSVC GCC
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) || defined(__WIN32__)
#ifndef SUANPAN_WIN
#define SUANPAN_WIN
#endif
#endif

// SUANPAN_WIN
#if defined(WIN64) || defined(_WIN64) || defined(__WIN64) || defined(__WIN64__)
#ifndef SUANPAN_WIN
#define SUANPAN_WIN
#endif
#endif

#ifdef SUANPAN_WIN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

// SUANPAN_UNIX
#if defined(unix) || defined(__unix__) || defined(__linux__) || defined(linux)
#ifndef SUANPAN_UNIX
#define SUANPAN_UNIX
#endif
#endif

#ifdef SUANPAN_VERSION
#undef SUANPAN_VERSION
#endif
#ifdef SUANPAN_COMPILER
#undef SUANPAN_COMPILER
#endif

// SUANPAN_VERSION SUANPAN_COMPILER
#ifdef __clang__
// __clang__ --> clang
#ifdef SUANPAN_VERSION
#undef SUANPAN_VERSION
#endif
#define SUANPAN_VERSION __VERSION__
#ifdef SUANPAN_COMPILER
#undef SUANPAN_COMPILER
#endif
#define SUANPAN_COMPILER "CLANG"
#define SUANPAN_CLANG
#elif defined(__GNUG__)
// __GNUG__ --> GCC
#define SUANPAN_VERSION __VERSION__
#define SUANPAN_COMPILER "GCC"
#define SUANPAN_GCC
#elif defined(_MSC_BUILD)
// _MSC_BUILD --> MSVC
#define SUANPAN_VERSION _MSC_FULL_VER
#define SUANPAN_COMPILER "MSVC"
#define SUANPAN_MSVC
// cuda unused local function
#pragma warning(disable : 4505)
#elif defined(__ICC)
// __ICC --> Intel C++
#define SUANPAN_VERSION __ICC
#define SUANPAN_COMPILER "INTEL"
#define SUANPAN_INTEL
#ifdef SUANPAN_WIN
#undef SUANPAN_WIN
#endif
#ifndef SUANPAN_UNIX
#define SUANPAN_UNIX
#endif
#elif defined(__ICL)
// __ICL --> Intel C++
#define SUANPAN_VERSION __ICL
#define SUANPAN_COMPILER "INTEL"
#define SUANPAN_INTEL
#ifdef SUANPAN_UNIX
#undef SUANPAN_UNIX
#endif
#ifndef SUANPAN_WIN
#define SUANPAN_WIN
#endif
#endif

// _USRDLL --> MSVC
#ifdef _USRDLL
#ifndef SUANPAN_DLL
#define SUANPAN_DLL
#endif
#endif

#ifdef SUANPAN_WIN
// WIN MSVC GCC IMPORT
#define SUANPAN_IMPORT extern "C" __declspec(dllimport)
// WIN MSVC GCC EXPORT
#define SUANPAN_EXPORT extern "C" __declspec(dllexport)
#elif defined(SUANPAN_UNIX)
// UNIX GCC IMPORT
#define SUANPAN_IMPORT extern "C"
// UNIX GCC EXPORT
#define SUANPAN_EXPORT extern "C"
#else
// EMPTY
#define SUANPAN_IMPORT extern "C"
#define SUANPAN_EXPORT extern "C"
#endif

#ifdef SUANPAN_DLL
SUANPAN_IMPORT bool SUANPAN_PRINT;
SUANPAN_IMPORT bool SUANPAN_VERBOSE;
SUANPAN_IMPORT bool SUANPAN_COLOR;
#else
SUANPAN_EXPORT bool SUANPAN_PRINT;
SUANPAN_EXPORT bool SUANPAN_VERBOSE;
SUANPAN_EXPORT bool SUANPAN_COLOR;
#endif

constexpr auto SUANPAN_EXIT = 1;
constexpr auto SUANPAN_SUCCESS = 0;
constexpr auto SUANPAN_FAIL = -1;

// TWO IMPLEMENTATIONS
#ifndef SUANPAN_WIN
#define _strcmpi strcasecmp
#endif

#ifdef SUANPAN_MT
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_sort.h>
#define suanpan_sort tbb::parallel_sort
#define suanpan_for_each tbb::parallel_for_each
#else
#define suanpan_sort std::sort
#define suanpan_for_each std::for_each
#endif

#include <iostream>
inline auto& SUANPAN_COUT = std::cout;
inline auto& SUANPAN_CWRN = std::cout;
inline auto& SUANPAN_CERR = std::cout;
inline auto& SUANPAN_CFTL = std::cout;

#define ARMA_COUT_STREAM SUANPAN_COUT
#define ARMA_CERR_STREAM SUANPAN_COUT

#include <armadillo/armadillo>
using namespace arma;

#include <filesystem>
namespace fs = std::filesystem;

#include <fmt/color.h>
#include <mutex>

namespace suanpan {
    inline std::mutex print_mutex;

    inline std::string pattern(const std::string_view header, const std::string_view& file_name, const std::string_view& format) {
        std::string pattern{header};
        pattern += fs::path(file_name).filename().string();
        pattern += ":{:<6d}";
        pattern += format;
        return pattern;
    }

    template<typename... T> void debug(const std::string_view file_name, const int line, const std::string_view format_str, const T&... args) {
        if(!SUANPAN_VERBOSE || !SUANPAN_PRINT) return;
        const std::scoped_lock lock(print_mutex);
        if(SUANPAN_COLOR) SUANPAN_COUT << fmt::vformat(fg(fmt::color::coral), pattern("[DEBUG] ", file_name, format_str), fmt::make_format_args(line, args...));
        else SUANPAN_COUT << fmt::vformat(pattern("[DEBUG] ", file_name, format_str), fmt::make_format_args(line, args...));
    }

    template<typename... T> void warning(const std::string_view file_name, const int line, const std::string_view format_str, const T&... args) {
        if(!SUANPAN_PRINT) return;
        const std::scoped_lock lock(print_mutex);
        if(SUANPAN_COLOR) SUANPAN_CWRN << fmt::vformat(fg(fmt::color::slate_blue), pattern("[WARNING] ", file_name, format_str), fmt::make_format_args(line, args...));
        else SUANPAN_CWRN << fmt::vformat(pattern("[WARNING] ", file_name, format_str), fmt::make_format_args(line, args...));
    }

    template<typename... T> void error(const std::string_view file_name, const int line, const std::string_view format_str, const T&... args) {
        if(!SUANPAN_PRINT) return;
        const std::scoped_lock lock(print_mutex);
        if(SUANPAN_COLOR) SUANPAN_CERR << fmt::vformat(fg(fmt::color::orange), pattern("[ERROR] ", file_name, format_str), fmt::make_format_args(line, args...));
        else SUANPAN_CERR << fmt::vformat(pattern("[ERROR] ", file_name, format_str), fmt::make_format_args(line, args...));
    }

    template<typename... T> void fatal(const std::string_view file_name, const int line, const std::string_view format_str, const T&... args) {
        if(!SUANPAN_PRINT) return;
        const std::scoped_lock lock(print_mutex);
        if(SUANPAN_COLOR) SUANPAN_CFTL << fmt::vformat(fg(fmt::color::violet), pattern("[FATAL] ", file_name, format_str), fmt::make_format_args(line, args...));
        else SUANPAN_CFTL << fmt::vformat(pattern("[FATAL] ", file_name, format_str), fmt::make_format_args(line, args...));
    }

    template<typename... T> void info(const std::string_view format_str, const T&... args) {
        if(!SUANPAN_PRINT) return;
        const std::scoped_lock lock(print_mutex);
        if(SUANPAN_COLOR) SUANPAN_COUT << fmt::vformat(fg(fmt::color::green_yellow), format_str, fmt::make_format_args(args...));
        else SUANPAN_COUT << fmt::vformat(format_str, fmt::make_format_args(args...));
    }

    template<typename... T> std::string format(const std::string_view format_str, const T&... args) { return fmt::vformat(format_str, fmt::make_format_args(args...)); }

    template<typename T> std::string format(const Col<T>& in_vec) {
        std::string output;
        if(std::is_floating_point_v<T>) for(const auto I : in_vec) output += format(" {: 1.4e}", I);
        else for(const auto I : in_vec) output += format(" {:6d}", I);
        output += '\n';
        return output;
    }

    template<typename T> void info(const Col<T>& in_vec) {
        if(!SUANPAN_PRINT) return;
        const std::scoped_lock lock(print_mutex);
        if(SUANPAN_COLOR) SUANPAN_COUT << fmt::format(fg(fmt::color::green_yellow), format(in_vec));
        else SUANPAN_COUT << format(in_vec);
    }

    template<typename T> void info(const std::string_view format_str, const Col<T>& in_vec) {
        if(!SUANPAN_PRINT) return;
        std::string output = format(format_str);
        if(format_str.back() != '\t' && format_str.back() != '\n') output += '\n';
        output += format(in_vec);
        const std::scoped_lock lock(print_mutex);
        if(SUANPAN_COLOR) SUANPAN_COUT << fmt::format(fg(fmt::color::green_yellow), output);
        else SUANPAN_COUT << output;
    }

    template<typename... T> void highlight(const std::string_view format_str, const T&... args) {
        if(!SUANPAN_PRINT) return;
        const std::scoped_lock lock(print_mutex);
        if(SUANPAN_COLOR) SUANPAN_COUT << fmt::vformat(fg(fmt::color::crimson), format_str, fmt::make_format_args(args...));
        else SUANPAN_COUT << fmt::vformat(format_str, fmt::make_format_args(args...));
    }
}

#ifdef SUANPAN_MSVC
#pragma warning(disable : 4100)
#endif
#include <functional>

inline void suanpan_assert(const std::function<void()>& F) {
#ifdef SUANPAN_DEBUG
    F();
#endif
}
#ifdef SUANPAN_MSVC
#pragma warning(default : 4100)
#endif

#define suanpan_info suanpan::info
#define suanpan_highlight suanpan::highlight
#define suanpan_debug(...) suanpan::debug(__FILE__, __LINE__, ##__VA_ARGS__)
#define suanpan_warning(...) suanpan::warning(__FILE__, __LINE__, ##__VA_ARGS__)
#define suanpan_error(...) suanpan::error(__FILE__, __LINE__, ##__VA_ARGS__)
#define suanpan_fatal(...) suanpan::fatal(__FILE__, __LINE__, ##__VA_ARGS__)

#include <memory>

using std::shared_ptr;
using std::unique_ptr;
using std::weak_ptr;

using std::make_shared;
using std::make_unique;

using std::exception;
using std::invalid_argument;
using std::logic_error;
using std::out_of_range;

using std::istringstream;
using std::ostringstream;
using std::string;

template<class T> concept sp_d = std::is_floating_point_v<T>;
template<class T> concept sp_i = std::is_integral_v<T>;

namespace suanpan {
    template<class IN, class FN> void for_all(IN& from, FN&& func) {
        suanpan_for_each(from.begin(), from.end(), std::forward<FN>(func));
    }
}

#if defined(SUANPAN_CLANG) && !defined(__cpp_lib_ranges)
// as of clang 13, ranges support is not complete
namespace std::ranges {
    template<class IN, class OUT, class FN> OUT transform(IN& from, OUT to, FN&& func) { return std::transform(from.begin(), from.end(), to, std::forward<FN>(func)); }

    template<class IN, class FN> FN for_each(IN& from, FN&& func) { return std::for_each(from.begin(), from.end(), std::forward<FN>(func)); }

    template<class IN, class OUT> OUT copy(IN& from, OUT to) { return std::copy(from.begin(), from.end(), to); }
} // namespace std::ranges
#endif

#endif
