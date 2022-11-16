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

#include <suanPan.h>
#include <Toolbox/argumentParser.h>

#ifdef SUANPAN_WIN
#include <Windows.h>

BOOL WIN_EVENT(DWORD) { return TRUE; }
#endif

// ReSharper disable once CppParameterMayBeConst
int main(int argc, char** argv) {
#ifdef SUANPAN_WIN
#if defined(SUANPAN_DEBUG) && defined(SUANPAN_MSVC)
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
    if(!SetConsoleCtrlHandler(WIN_EVENT, TRUE)) return 0;
    const auto handle = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO info;
    GetConsoleScreenBufferInfo(handle, &info);
    const auto current_attribute = info.wAttributes;
    SetConsoleTextAttribute(handle, FOREGROUND_INTENSITY | FOREGROUND_GREEN);
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#else
    SUANPAN_SYNC_COUT << FOREGROUND_GREEN;
#endif

#ifdef SUANPAN_DEBUG
    argument_parser(argc, argv);
#else
    try { argument_parser(argc, argv); }
    catch(const std::exception& e) { suanpan_fatal("some unexpected error happens: %s, please file a bug report via https://github.com/TLCFEM/suanPan/issues.\n", e.what()); }
#endif

#ifdef SUANPAN_WIN
    SetConsoleTextAttribute(handle, current_attribute);
#else
    SUANPAN_SYNC_COUT << "\033[0m";
#endif

    return SUANPAN_SUCCESS;
}
