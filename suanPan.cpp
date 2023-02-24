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

#include <Toolbox/argument.h>
#include <suanPan.h>

#ifdef SUANPAN_WIN
#include <Windows.h>
#endif

#ifdef SUANPAN_MAGMA
#include <magma_auxiliary.h>
#endif

// ReSharper disable once CppParameterMayBeConst
int main(int argc, char** argv) {
#ifdef SUANPAN_WIN
#if defined(SUANPAN_DEBUG) && defined(SUANPAN_MSVC)
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif

#ifdef SUANPAN_MAGMA
    magma_init();
#endif

#ifdef SUANPAN_DEBUG
#ifdef SUANPAN_MAGMA
    magma_print_environment();
#endif
    argument_parser(argc, argv);
#else
    try { argument_parser(argc, argv); }
    catch(const std::exception& e) { suanpan_fatal("Some unexpected error happens: {}, please file a bug report via https://github.com/TLCFEM/suanPan/issues.\n", e.what()); }
#endif

#ifdef SUANPAN_MAGMA
    magma_finalize();
#endif

    return SUANPAN_SUCCESS;
}
