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

#include "argumentParser.h"
#include <Step/Bead.h>
#include <Toolbox/Converter.h>
#include <Toolbox/commandParser.h>
#include <Toolbox/utility.h>
#include <UnitTest/CatchTest.h>
#include <iomanip>
#include "revision.h"
#ifdef SUANPAN_WIN
#include <Windows.h>
#endif

#ifdef SUANPAN_MSVC
#pragma warning(disable : 4702 4996)
#endif

using std::ifstream;
using std::ofstream;

constexpr auto SUANPAN_MAJOR = 2;
constexpr auto SUANPAN_MINOR = 8;
constexpr auto SUANPAN_PATCH = 0;
constexpr auto SUANPAN_CODE = "Betelgeuse";
constexpr auto SUANPAN_ARCH = 64;

bool SUANPAN_PRINT = true;
#ifdef SUANPAN_DEBUG
bool SUANPAN_VERBOSE = true;
#else
bool SUANPAN_VERBOSE = false;
#endif

const char* SUANPAN_EXE = "";

bool check_debugger() {
#ifdef SUANPAN_DEBUG
    return false;
#endif
#ifdef SUANPAN_WIN
    if(IsDebuggerPresent()) exit(EXIT_SUCCESS); // NOLINT(concurrency-mt-unsafe)

    BOOL FLAG = false;
    if(CheckRemoteDebuggerPresent(GetCurrentProcess(), &FLAG) && FLAG) std::swap(access::rw(datum::eps), access::rw(datum::pi));

    CONTEXT CTX{};
    CTX.ContextFlags = CONTEXT_DEBUG_REGISTERS;
    if(GetThreadContext(GetCurrentThread(), &CTX) && (CTX.Dr0 != 0 || CTX.Dr1 != 0 || CTX.Dr2 != 0 || CTX.Dr3 != 0)) exit(EXIT_SUCCESS); // NOLINT(concurrency-mt-unsafe)
#elif defined(SUANPAN_UNIX)
#endif
    return false;
}

void check_version(const char* path_to_executable) {
    auto updater_module = fs::path(path_to_executable).parent_path();

#ifdef SUANPAN_WIN
    updater_module.append("updater.exe");
#else
    updater_module.append("updater");
#endif

    if(!exists(updater_module)) return;

    auto terminal = istringstream("\"" + updater_module.string() + "\" " + std::to_string(100 * SUANPAN_MAJOR + 10 * SUANPAN_MINOR + SUANPAN_PATCH));

    execute_command(terminal);
}

void strip_mode(const string& input_file_name, const string& output_file_name) {
    ifstream input_file(input_file_name);

    if(!input_file.is_open()) {
        suanpan_error("Fail to open \"{}\".\n", input_file_name);
        return;
    }

    ofstream output_file(output_file_name);

    if(!output_file.is_open()) {
        suanpan_error("Fail to open \"{}\".\n", output_file_name);
        return;
    }

    output_file.setf(std::ios::scientific);
    output_file << std::setprecision(3);

    string line;

    while(std::getline(input_file, line)) {
        if(line.empty() || if_contain(line, "**")) continue;

        for(auto& I : line) I = static_cast<char>(std::tolower(static_cast<int>(I)));

        output_file << line << '\n';
    }
}

void convert_mode(const string& input_file_name, const string& output_file_name) {
    ifstream input_file(input_file_name);

    if(!input_file.is_open()) {
        suanpan_error("Fail to open \"{}\".\n", input_file_name);
        return;
    }

    ofstream output_file(output_file_name);

    if(!output_file.is_open()) {
        suanpan_error("Fail to open \"{}\".\n", output_file_name);
        return;
    }

    const auto pos = output_file_name.find_last_of('/');

    Converter abaqus_converter(string::npos == pos ? "" : output_file_name.substr(0, pos + 1));

    abaqus_converter.process(input_file, output_file);
}

void print_header() {
    suanpan_info("+--------------------------------------------------+\n");
    suanpan_info("|   __        __         suanPan is an open source |\n");
    suanpan_info("|  /  \\      |  \\           FEM framework ({}-bit) |\n", SUANPAN_ARCH);
    suanpan_info("|  \\__       |__/  __   __      {} ({}.{}.{}) |\n", SUANPAN_CODE, SUANPAN_MAJOR, SUANPAN_MINOR, SUANPAN_PATCH);
    suanpan_info("|     \\ |  | |    |  \\ |  |      by tlc @ {} |\n", SUANPAN_REVISION);
    suanpan_info("|  \\__/ |__| |    |__X |  |    all rights reserved |\n");
    suanpan_info("|                           10.5281/zenodo.1285221 |\n");
    suanpan_info("+--------------------------------------------------+\n");
#ifdef SUANPAN_WIN
    suanpan_info("|  https://github.com/TLCFEM/suanPan               |\n");
    suanpan_info("|  https://github.com/TLCFEM/suanPan-manual        |\n");
    suanpan_info("+--------------------------------------------------+\n");
    suanpan_info("|  https://gitter.im/suanPan-dev/community         |\n");
#else
    static constexpr std::array POOL{"\U0001F308", "\U0001F30F", "\U0001F3A7", "\U0001F3B1", "\U0001F479", "\U0001F4BB", "\U0001F50B", "\U0001F514", "\U0001F680", "\U0001F9E9"};
    arma_rng::set_seed_random();
    suanpan_info("|  \U0001F9EE https://github.com/TLCFEM/suanPan            |\n");
    suanpan_info("|  \U0001F4DA https://github.com/TLCFEM/suanPan-manual     |\n");
    suanpan_info("+--------------------------------------------------+\n");
    suanpan_info("|  {} https://gitter.im/suanPan-dev/community      |\n", POOL[randi() % POOL.size()]);
#endif
    suanpan_info("+--------------------------------------------------+\n\n");
}

void argument_parser(const int argc, char** argv) {
    if(check_debugger()) return;

    if(0 == argc) return;

    SUANPAN_EXE = argv[0];

    string input_file_name;
    const auto buffer_backup = SUANPAN_COUT.rdbuf();

    wall_clock T;
    T.tic();

    if(argc > 1) {
        ofstream output_file;
        string output_file_name;
        auto strip = false, convert = false, check_new = true;

        for(auto I = 1; I < argc; ++I) {
            if(is_equal(argv[I], "-v") || is_equal(argv[I], "--version")) print_version();
            else if(is_equal(argv[I], "-h") || is_equal(argv[I], "--help")) print_helper();
            else if(is_equal(argv[I], "-t") || is_equal(argv[I], "--test")) test_mode();
            else if(is_equal(argv[I], "-f") || is_equal(argv[I], "--file")) input_file_name = argv[++I];
            else if(is_equal(argv[I], "-o") || is_equal(argv[I], "--output")) output_file_name = argv[++I];
            else if(is_equal(argv[I], "-np") || is_equal(argv[I], "--noprint")) SUANPAN_PRINT = false;
            else if(is_equal(argv[I], "-vb") || is_equal(argv[I], "--verbose")) SUANPAN_VERBOSE = true;
            else if(is_equal(argv[I], "-nu") || is_equal(argv[I], "--noupdate")) check_new = false;
            else if(is_equal(argv[I], "-s") || is_equal(argv[I], "--strip")) {
                strip = true;
                convert = false;
            }
            else if(is_equal(argv[I], "-c") || is_equal(argv[I], "--convert")) {
                convert = true;
                strip = false;
            }
            else if(is_equal(argv[I], "-ctest") || is_equal(argv[I], "--catch2test")) {
                catchtest_main(argc, argv);
                return;
            }
        }

        if(check_new) check_version(SUANPAN_EXE);

        if(strip || convert) {
            if(input_file_name.empty()) return;

            if(output_file_name.empty()) {
                output_file_name = input_file_name;

                auto found = output_file_name.rfind(".inp");

                if(string::npos == found) found = output_file_name.rfind(".INP");

                if(string::npos != found) output_file_name.erase(found, 4);

                output_file_name += strip ? "_out.inp" : "_out.supan";
            }

            for(auto& I : output_file_name) if(I == '\\') I = '/';

            return convert ? convert_mode(input_file_name, output_file_name) : strip_mode(input_file_name, output_file_name);
        }

        if(!output_file_name.empty()) {
            output_file.open(output_file_name);
            if(output_file.is_open()) SUANPAN_COUT.rdbuf(output_file.rdbuf());
            else
                suanpan_error("Cannot open the output file \"{}\".\n", output_file_name);
        }

        if(!input_file_name.empty()) {
            print_header();
            if(const auto model = make_shared<Bead>(); process_file(model, input_file_name.c_str()) != SUANPAN_EXIT) {
                if(output_file.is_open()) {
                    SUANPAN_COUT.rdbuf(buffer_backup);
                    print_header();
                }
                cli_mode(model);
            }
        }
    }
    else {
        check_version(SUANPAN_EXE);
        print_header();
        const auto model = make_shared<Bead>();
        cli_mode(model);
    }

    suanpan_info("\nTime Wasted: {:.4f} Seconds.\n", T.toc());
}

void print_version() {
    suanpan_info("Copyright (C) 2017-2023 Theodore Chang\n\n");
    suanpan_info("This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\n\n");
    suanpan_info("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.\n\n");
    suanpan_info("You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\n\n");
    suanpan_info("suanPan is an open source FEM framework.\n");
    suanpan_info("    The binary is compiled on {}.\n", __DATE__);
    suanpan_info("    The source code of suanPan is hosted on GitHub. https://tlcfem.github.io/suanPan/\n");
    suanpan_info("    The documentation is hosted on GitBook and readthedocs. https://tlcfem.gitbook.io/suanpan-manual/ and https://suanpan-manual.readthedocs.io/\n");
#ifdef SUANPAN_MKL
    suanpan_info("    The linear algebra support is provided by Armadillo with Intel MKL. http://arma.sourceforge.net/\n");
#else
    suanpan_info("    The linear algebra support is provided by Armadillo. http://arma.sourceforge.net/\n");
#endif
#ifdef SUANPAN_CUDA
    suanpan_info("    The GPCPU solvers are provided by CUDA. https://developer.nvidia.com/about-cuda\n");
#endif
#ifdef SUANPAN_MT
    suanpan_info("    The parallelisation support is implemented via TBB library. https://github.com/oneapi-src/oneTBB\n");
#endif
#ifdef SUANPAN_VTK
    suanpan_info("    The visualisation support is implemented via VTK library. https://vtk.org/\n");
#endif
    suanpan_info("\nPlease join gitter for any feedback. https://gitter.im/suanPan-dev/community\n");
    suanpan_info("\n\n[From Wikipedia] Betelgeuse is usually the tenth-brightest star in the night sky and, after Rigel, the second-brightest in the constellation of Orion. It is a distinctly reddish semiregular variable star whose apparent magnitude has the widest range displayed by any first-magnitude star.\n\n");
}

void print_helper() {
    suanpan_info("Available Parameters:\n");
    suanpan_info("\t-{:<10}  --{:<20}{}\n", "v", "version", "check version information");
    suanpan_info("\t-{:<10}  --{:<20}{}\n", "h", "help", "print this helper");
    suanpan_info("\t-{:<10}  --{:<20}{}\n", "s", "strip", "strip comments out in given ABAQUS input file");
    // suanpan_info("\t-{:<10}  --{:<20}{}\n", "c", "convert", "partially convert ABAQUS input file into suanPan model script");
    suanpan_info("\t-{:<10}  --{:<20}{}\n", "np", "noprint", "suppress most console output");
    suanpan_info("\t-{:<10}  --{:<20}{}\n", "nu", "noupdate", "do not check for newer version on startup");
    suanpan_info("\t-{:<10}  --{:<20}{}\n", "f", "file", "process model file");
    suanpan_info("\t-{:<10}  --{:<20}{}\n", "o", "output", "set output file for logging");
    suanpan_info("\n");
}

void cli_mode(const shared_ptr<Bead>& model) {
#ifdef SUANPAN_WIN
    // ReSharper disable once CppDeprecatedEntity
    auto history_path = fs::path(getenv("USERPROFILE")); // NOLINT(concurrency-mt-unsafe, clang-diagnostic-deprecated-declarations)
#else
    auto history_path = fs::path(getenv("HOME"));
#endif

    history_path.append(".suanpan-history.sp");

    ofstream output_file(history_path, std::ios_base::app | std::ios_base::out);

    string all_line;
    while(true) {
        string command_line;
        suanpan_info("suanPan ~<> ");
        getline(std::cin, command_line);
        if(!command_line.empty() && command_line[0] != '#' && command_line[0] != '!') {
            if(const auto if_comment = command_line.find('!'); string::npos != if_comment) command_line.erase(if_comment);
            for(auto& c : command_line) if(',' == c || '\t' == c || '\r' == c || '\n' == c) c = ' ';
            while(!command_line.empty() && *command_line.crbegin() == ' ') command_line.pop_back();
            if(command_line.empty()) continue;
            if(*command_line.crbegin() == '\\') {
                command_line.back() = ' ';
                all_line.append(command_line);
            }
            else {
                all_line.append(command_line);
                istringstream tmp_str(all_line);
                if(output_file.is_open()) output_file << all_line << '\n';
                if(process_command(model, tmp_str) == SUANPAN_EXIT) return;
                all_line.clear();
            }
        }
    }
}
