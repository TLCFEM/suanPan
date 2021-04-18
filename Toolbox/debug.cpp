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

#include "debug.h"
#include <cstdarg>
#include <iostream>
#include <suanPan.h>
#include <vector>
#ifdef SUANPAN_WIN
#include <Windows.h>
#endif

#ifdef SUANPAN_MSVC
#pragma warning(disable:4702 4100)
#endif

void suanpan_info(const char* M, ...) {
	if(!SUANPAN_PRINT) return;

	va_list arguments_a, arguments_b;
	va_start(arguments_a, M);
	va_copy(arguments_b, arguments_a);
	std::vector<char> buffer(1ll + vsnprintf(nullptr, 0, M, arguments_a));
	va_end(arguments_a);
	vsnprintf(&buffer[0], buffer.size(), M, arguments_b);
	va_end(arguments_b);
	std::cout << buffer.data();
}

// empty function call will automatically be optimized out
void suanpan_debug(const char* M, ...) {
#ifdef SUANPAN_DEBUG
#ifdef SUANPAN_WIN
	const auto handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(handle, &info);
	const auto current_attribute = info.wAttributes;
	SetConsoleTextAttribute(handle, FOREGROUND_INTENSITY | FOREGROUND_GREEN);
#else
	std::cout << "\033[1;36m";
#endif
	cout << "debug: ";
#ifdef SUANPAN_WIN
	SetConsoleTextAttribute(handle, current_attribute);
#else
	std::cout << "\033[1;32m";
#endif
	va_list arguments_a, arguments_b;
	va_start(arguments_a, M);
	va_copy(arguments_b, arguments_a);
	std::vector<char> buffer(1ll + vsnprintf(nullptr, 0, M, arguments_a));
	va_end(arguments_a);
	vsnprintf(&buffer[0], buffer.size(), M, arguments_b);
	va_end(arguments_b);
	std::cout << buffer.data();
#endif
}

// empty function call will automatically be optimized out
void suanpan_extra_debug(const char* M, ...) {
#ifdef SUANPAN_EXTRA_DEBUG
#ifdef SUANPAN_WIN
	const auto handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(handle, &info);
	const auto current_attribute = info.wAttributes;
	SetConsoleTextAttribute(handle, FOREGROUND_INTENSITY | FOREGROUND_GREEN);
#else
	std::cout << "\033[1;36m";
#endif
	cout << "extra debug: ";
#ifdef SUANPAN_WIN
	SetConsoleTextAttribute(handle, current_attribute);
#else
	std::cout << "\033[1;32m";
#endif
	va_list arguments_a, arguments_b;
	va_start(arguments_a, M);
	va_copy(arguments_b, arguments_a);
	std::vector<char> buffer(1ll + vsnprintf(nullptr, 0, M, arguments_a));
	va_end(arguments_a);
	vsnprintf(&buffer[0], buffer.size(), M, arguments_b);
	va_end(arguments_b);
#ifndef SUANPAN_WIN
	std::cout << "\033[1;32m";
#endif
	std::cout << buffer.data();
#endif
}

void suanpan_warning(const char* M, ...) {
#ifdef SUANPAN_WIN
	const auto handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(handle, &info);
	const auto current_attribute = info.wAttributes;
	SetConsoleTextAttribute(handle, FOREGROUND_INTENSITY | FOREGROUND_BLUE);
#else
	std::cout << "\033[1;34m";
#endif
	std::cout << "warning: ";
#ifdef SUANPAN_WIN
	SetConsoleTextAttribute(handle, current_attribute);
#else
	std::cout << "\033[1;32m";
#endif
	va_list arguments_a, arguments_b;
	va_start(arguments_a, M);
	va_copy(arguments_b, arguments_a);
	std::vector<char> buffer(1ll + vsnprintf(nullptr, 0, M, arguments_a));
	va_end(arguments_a);
	vsnprintf(&buffer[0], buffer.size(), M, arguments_b);
	va_end(arguments_b);
	std::cout << buffer.data();
}

void suanpan_error(const char* M, ...) {
#ifdef SUANPAN_WIN
	const auto handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(handle, &info);
	const auto current_attribute = info.wAttributes;
	SetConsoleTextAttribute(handle, FOREGROUND_INTENSITY | FOREGROUND_GREEN | FOREGROUND_RED);
#else
	std::cout << "\033[1;33m";
#endif
	std::cout << "error: ";
#ifdef SUANPAN_WIN
	SetConsoleTextAttribute(handle, current_attribute);
#else
	std::cout << "\033[1;32m";
#endif
	va_list arguments_a, arguments_b;
	va_start(arguments_a, M);
	va_copy(arguments_b, arguments_a);
	std::vector<char> buffer(1ll + vsnprintf(nullptr, 0, M, arguments_a));
	va_end(arguments_a);
	vsnprintf(&buffer[0], buffer.size(), M, arguments_b);
	va_end(arguments_b);
	std::cout << buffer.data();
}

void suanpan_fatal(const char* M, ...) {
#ifdef SUANPAN_WIN
	const auto handle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(handle, &info);
	const auto current_attribute = info.wAttributes;
	SetConsoleTextAttribute(handle, FOREGROUND_INTENSITY | FOREGROUND_RED);
#else
	std::cout << "\033[1;31m";
#endif
	std::cout << "fatal: ";
#ifdef SUANPAN_WIN
	SetConsoleTextAttribute(handle, current_attribute);
#else
	std::cout << "\033[1;32m";
#endif
	va_list arguments_a, arguments_b;
	va_start(arguments_a, M);
	va_copy(arguments_b, arguments_a);
	std::vector<char> buffer(1ll + vsnprintf(nullptr, 0, M, arguments_a));
	va_end(arguments_a);
	vsnprintf(&buffer[0], buffer.size(), M, arguments_b);
	va_end(arguments_b);
	std::cout << buffer.data();
}

void suanpan_print(const char* M, ...) {
	va_list arguments_a, arguments_b;
	va_start(arguments_a, M);
	va_copy(arguments_b, arguments_a);
	std::vector<char> buffer(1ll + vsnprintf(nullptr, 0, M, arguments_a));
	va_end(arguments_a);
	vsnprintf(&buffer[0], buffer.size(), M, arguments_b);
	va_end(arguments_b);
	std::cout << buffer.data();
}

void suanpan_debug(const std::function<void()>& F) {
#ifdef SUANPAN_DEBUG
	F();
#endif
}

bool check_debugger() {
#ifdef SUANPAN_DEBUG
	return false;
#endif
#ifdef SUANPAN_WIN
	if(IsDebuggerPresent()) exit(EXIT_SUCCESS);

	BOOL FLAG = false;
	if(CheckRemoteDebuggerPresent(GetCurrentProcess(), &FLAG) && FLAG) swap_constant();

	CONTEXT CTX{};
	CTX.ContextFlags = CONTEXT_DEBUG_REGISTERS;
	if(GetThreadContext(GetCurrentThread(), &CTX) && (CTX.Dr0 != 0 || CTX.Dr1 != 0 || CTX.Dr2 != 0 || CTX.Dr3 != 0)) exit(EXIT_SUCCESS);

#elif defined(SUANPAN_UNIX)
#endif
	return false;
}

void swap_constant() { std::swap(access::rw(datum::eps), access::rw(datum::pi)); }

#ifdef SUANPAN_MSVC
#pragma warning(default:4702 4100)
#endif
