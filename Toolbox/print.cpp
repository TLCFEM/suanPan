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

#include "print.h"
#include <cstdarg>
#include <suanPan.h>
#ifdef SUANPAN_WIN
#include <Windows.h>
#endif

#ifdef SUANPAN_MSVC
#pragma warning(disable : 4100)
#endif

#ifdef SUANPAN_WIN
void suanpan_print_header(const char* header, const unsigned short color) {
    const auto handle = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO info;
    GetConsoleScreenBufferInfo(handle, &info);
    const auto current_attribute = info.wAttributes;
    SetConsoleTextAttribute(handle, FOREGROUND_INTENSITY | color);
    SUANPAN_SYNC_COUT << header << ": ";
    SetConsoleTextAttribute(handle, current_attribute);
}
#else
void suanpan_print_header(const char* header, const char* color) { SUANPAN_SYNC_COUT << color << header << ": " << FOREGROUND_GREEN; }
#endif

namespace suanpan {
    std::mutex print_mutex;
}

void suanpan_info(const char* M, ...) {
    if(!SUANPAN_PRINT) return;

    va_list list;
    va_start(list, M);
    const size_t len = std::vsnprintf(nullptr, 0, M, list);
    va_end(list);
    std::vector<char> vec(len + 1);
    va_start(list, M);
    std::vsnprintf(vec.data(), len + 1, M, list);
    va_end(list);

    const std::scoped_lock lock(suanpan::print_mutex);
    SUANPAN_SYNC_COUT << vec.data();
}

// empty function call will automatically be optimized out
void suanpan_debug(const char* M, ...) {
#ifdef SUANPAN_DEBUG
    va_list list;
    va_start(list, M);
    const size_t len = std::vsnprintf(nullptr, 0, M, list);
    va_end(list);
    std::vector<char> vec(len + 1);
    va_start(list, M);
    std::vsnprintf(vec.data(), len + 1, M, list);
    va_end(list);

    const std::scoped_lock lock(suanpan::print_mutex);
    suanpan_print_header("debug", FOREGROUND_CYAN);
    SUANPAN_SYNC_COUT << vec.data();
#endif
}

// empty function call will automatically be optimized out
void suanpan_extra_debug(const char* M, ...) {
#ifdef SUANPAN_EXTRA_DEBUG
    if(!SUANPAN_VERBOSE) return;

    va_list list;
    va_start(list, M);
    const size_t len = std::vsnprintf(nullptr, 0, M, list);
    va_end(list);
    std::vector<char> vec(len + 1);
    va_start(list, M);
    std::vsnprintf(vec.data(), len + 1, M, list);
    va_end(list);

    const std::scoped_lock lock(suanpan::print_mutex);
    suanpan_print_header("extra debug", FOREGROUND_CYAN);
    SUANPAN_SYNC_COUT << vec.data();
#endif
}

void suanpan_warning(const char* M, ...) {
    va_list list;
    va_start(list, M);
    const size_t len = std::vsnprintf(nullptr, 0, M, list);
    va_end(list);
    std::vector<char> vec(len + 1);
    va_start(list, M);
    std::vsnprintf(vec.data(), len + 1, M, list);
    va_end(list);

    const std::scoped_lock lock(suanpan::print_mutex);
    suanpan_print_header("warning", FOREGROUND_BLUE);
    SUANPAN_SYNC_COUT << vec.data();
}

void suanpan_error(const char* M, ...) {
    va_list list;
    va_start(list, M);
    const size_t len = std::vsnprintf(nullptr, 0, M, list);
    va_end(list);
    std::vector<char> vec(len + 1);
    va_start(list, M);
    std::vsnprintf(vec.data(), len + 1, M, list);
    va_end(list);

    const std::scoped_lock lock(suanpan::print_mutex);
    suanpan_print_header("error", FOREGROUND_YELLOW);
    SUANPAN_SYNC_COUT << vec.data();
}

void suanpan_fatal(const char* M, ...) {
    va_list list;
    va_start(list, M);
    const size_t len = std::vsnprintf(nullptr, 0, M, list);
    va_end(list);
    std::vector<char> vec(len + 1);
    va_start(list, M);
    std::vsnprintf(vec.data(), len + 1, M, list);
    va_end(list);

    const std::scoped_lock lock(suanpan::print_mutex);
    suanpan_print_header("fatal", FOREGROUND_RED);
    SUANPAN_SYNC_COUT << vec.data();
}

void suanpan_debug(const std::function<void()>& F) {
#ifdef SUANPAN_DEBUG
    F();
#endif
}
