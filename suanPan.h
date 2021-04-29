/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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

constexpr auto SUANPAN_MAJOR = 1;
constexpr auto SUANPAN_MINOR = 5;
constexpr auto SUANPAN_PATCH = 0;
constexpr auto SUANPAN_CODE = "     Acrux";
constexpr auto SUANPAN_ARCH = 64;

// SUANPAN_DEBUG
// _DEBUG --> MSVC
// DEBUG --> GCC
#if defined(_DEBUG) || defined(DEBUG)
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
// __GNUG__ --> GCC
#ifdef __GNUG__
#define SUANPAN_VERSION __VERSION__
#define SUANPAN_COMPILER "GCC"
#define SUANPAN_GCC
#endif
// _MSC_BUILD --> MSVC
#ifdef _MSC_BUILD
#define SUANPAN_VERSION _MSC_FULL_VER
#define SUANPAN_COMPILER "MSVC"
#define SUANPAN_MSVC
// cuda unused local function
#pragma warning(disable : 4505)
#endif
// __clang__ --> clang
#ifdef __clang__
#ifdef SUANPAN_VERSION
#undef SUANPAN_VERSION
#endif
#define SUANPAN_VERSION __VERSION__
#ifdef SUANPAN_COMPILER
#undef SUANPAN_COMPILER
#endif
#define SUANPAN_COMPILER "CLANG"
#define SUANPAN_CLANG
#endif
// __ICC --> Intel C++
#ifdef __ICC
#define SUANPAN_VERSION __ICC
#define SUANPAN_COMPILER "INTEL"
#define SUANPAN_INTEL
#ifdef SUANPAN_WIN
#undef SUANPAN_WIN
#endif
#ifndef SUANPAN_UNIX
#define SUANPAN_UNIX
#endif
#endif
// __ICL --> Intel C++
#ifdef __ICL
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
#else
SUANPAN_EXPORT bool SUANPAN_PRINT;
SUANPAN_EXPORT const char* SUANPAN_EXE;
#endif

#ifndef SUANPAN_EXIT
#define SUANPAN_EXIT 1
#endif

#ifndef SUANPAN_SUCCESS
#define SUANPAN_SUCCESS 0
#endif

#ifndef SUANPAN_FAIL
#define SUANPAN_FAIL (-1)
#endif

// TWO IMPLEMENTATIONS
#ifndef SUANPAN_WIN
#define _strcmpi strcasecmp
#endif

#include <Toolbox/debug.h>

#include <armadillo/armadillo>
#include <filesystem>
#include <future>
#include <memory>
#include <typeinfo>
namespace fs = std::filesystem;

#ifdef SUANPAN_MT
// ppl as alternative for MSVC platform
// #include <ppl.h>
// #define suanpan_for_each Concurrency::parallel_for_each
#ifdef SUANPAN_MSVC
#include <tbb/tbbmalloc_proxy.h>
#endif
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_sort.h>
#define suanpan_sort tbb::parallel_sort
#define suanpan_for_each tbb::parallel_for_each
#else
#define suanpan_sort std::sort
#define suanpan_for_each std::for_each
#endif

using namespace arma;

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

#endif
