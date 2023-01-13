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

// ReSharper disable StringLiteralTypo
// ReSharper disable IdentifierTypo
#include "command.h"
#include <Constraint/Constraint.h>
#include <Constraint/ConstraintParser.h>
#include <Converger/Converger.h>
#include <Converger/ConvergerParser.h>
#include <Domain/Domain.h>
#include <Domain/ExternalModule.h>
#include <Domain/Group/ElementGroup.h>
#include <Domain/Group/GroupGroup.h>
#include <Domain/Group/NodeGroup.h>
#include <Domain/Node.h>
#include <Element/Element.h>
#include <Element/ElementParser.h>
#include <Element/Visualisation/vtkParser.h>
#include <Load/Amplitude/Amplitude.h>
#include <Load/Load.h>
#include <Load/LoadParser.h>
#include <Material/Material.h>
#include <Material/MaterialParser.h>
#include <Material/MaterialTester.h>
#include <Recorder/Recorder.h>
#include <Recorder/RecorderParser.h>
#include <Section/SectionParser.h>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Solver.h>
#include <Solver/SolverParser.h>
#include <Step/Bead.h>
#include <Step/Frequency.h>
#include <Step/StepParser.h>
#include <Toolbox/Expression.h>
#include <Toolbox/ExpressionParser.h>
#include <Toolbox/argument.h>
#include <Toolbox/resampling.h>
#include <Toolbox/response_spectrum.h>
#include <Toolbox/thread_pool.hpp>
#include <thread>
#ifdef SUANPAN_MKL
#include <Domain/MetaMat/SparseMatFGMRES.hpp>
#endif
#ifdef SUANPAN_WIN
#include <Windows.h>
#endif

using std::ifstream;
using std::string;
using std::vector;

int SUANPAN_NUM_THREADS = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
fs::path SUANPAN_OUTPUT = fs::current_path();
extern fs::path SUANPAN_EXE;

void qrcode() {
    for(constexpr char encode[] = "SLLLLLLLWWWLWWWLWWWLWWWLLLLLLLSFWLLLWFWLUWLWUWLWWFFFWFWLLLWFSFWFFFWFWWFWWFFWWFUFUWWFWFFFWFSFLLLLLFWLWFUFWFUFUFULWFLLLLLFSLLLWLLLLFWWULWWULUUFFLLWWWLWWSULUUFFLWWULFFULFFWWUFLFWLULLFSLUUFWULFWUFLUUFLFFFUULLUULWFLSLUFULULLWUUUWLUULLWUUUFWLFWLFSLFLLLLLWLFWULWWLFFULFUFLWFWFLSLWLWWULLFWLFFULWUFFWWFULLUULFSLULFUFLFFFFLUUFULFUFFFFFFUWUWSLLLLLLLWFLUUWLUWFUUFFWLWFLUFFSFWLLLWFWFFWULWWUWFUWFLLLFUWWLSFWFFFWFWLFWFFULUFULLUWWFFLUUFSFLLLLLFWFFFLUUFLFFUFFFWLFWWFL"; const auto I : encode)
        if(I == 'S')
            suanpan_info("\n            ");
        else if(I == 'W')
            suanpan_info(" ");
        else if(I == 'F')
            suanpan_info("\xE2\x96\x88");
        else if(I == 'L')
            suanpan_info("\xE2\x96\x84");
        else if(I == 'U')
            suanpan_info("\xE2\x96\x80");

    suanpan_info("\n\n");
}

int benchmark() {
    constexpr auto N = 50;
    constexpr auto M = 5120;

    thread_pool pool(1);

    const mat A = mat(M, M, fill::randu) + eye(M, M);
    const vec b(M, fill::randu);

    const auto start = std::chrono::high_resolution_clock::now();

    for(auto I = 1; I <= N; ++I) {
        pool.push_task([I] {
            SUANPAN_COUT << '[';
            const auto length = static_cast<int>(50. * I / N);
            for(auto J = 0; J < length; ++J) SUANPAN_COUT << '=';
            for(auto J = length; J < 50; ++J) SUANPAN_COUT << '-';
            SUANPAN_COUT << "]\r";
            SUANPAN_COUT.flush();
        });
        vec x = solve(A, b);
        x(randi<uvec>(1, distr_param(0, M - 1))).fill(I);
    }

    const auto end = std::chrono::high_resolution_clock::now();

    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    pool.wait_for_tasks();

    suanpan_info("\nCurrent platform rates (higher is better): {:.2f}.\n", 1E9 / static_cast<double>(duration.count()));

    return SUANPAN_SUCCESS;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"
void overview() {
    const auto new_model = make_shared<Bead>();

    auto guide_command = [&](const string& target) {
        while(true) {
            string command_line;
            suanpan_highlight("overview -> ");
            getline(std::cin, command_line);
            if(is_equal(command_line, "q")) {
                suanpan_info("Bye. Happy modelling!\n\n");
                return SUANPAN_EXIT;
            }
            if(is_equal(command_line, target)) {
                istringstream tmp_str(command_line);
                const auto code = process_command(new_model, tmp_str);
                suanpan_highlight("overview -> [Press Enter to Continue]");
                getline(std::cin, command_line);
                return code;
            }
            suanpan_info("Try again.\n");
        }
    };

    std::streambuf* buffer_backup;
    unique_ptr<std::stringbuf> buffer_tmp;

    const auto redirect = [&] {
        buffer_backup = SUANPAN_COUT.rdbuf();
        buffer_tmp = make_unique<std::stringbuf>();
        SUANPAN_COUT.rdbuf(buffer_tmp.get());
    };
    const auto restore = [&] {
        SUANPAN_COUT.rdbuf(buffer_backup);
        for(const auto I : buffer_tmp->str()) {
            SUANPAN_COUT << I;
            SUANPAN_COUT.flush();
            std::this_thread::sleep_for(std::chrono::milliseconds(8));
        }
    };

    redirect();
    suanpan_info("Welcome to suanPan, a parallel finite element analysis software.\n\n"
        "It can be used to solve various types of solid mechanics problems. It also provides a flexible platform for researchers to develop new numerical models and algorithms in a hassle-free manner.\n\n"
        "The following is a quick introduction of the application, regarding its functionalities and basic usage. At any time, type in '");
    suanpan_highlight("q");
    suanpan_info("' to quit this introduction. First of all, type in '");
    suanpan_highlight("version");
    suanpan_info("' to check the license and version information.\n");
    restore();

    if(SUANPAN_EXIT == guide_command("version")) return;

    redirect();
    suanpan_info("By invoking the application without input file, you are currently in the interactive mode, which allows you to create finite element models interactively. "
        "A numerical model typically consists of nodes, elements connecting nodes, materials attached to elements, loads applied to nodes, boundary conditions imposed on nodes, etc. "
        "To analyze, a number of analysis steps need to be defined with the proper convergence criteria and solvers.\n\nType in '");
    suanpan_highlight("example");
    suanpan_info("' to check out a simple elastic cantilever example.\n");
    restore();

    if(SUANPAN_EXIT == guide_command("example")) return;

    redirect();
    suanpan_info("If you have some experience with ABAQUS scripting, you may have realised that the structure resembles that of ABAQUS *.inp files. "
        "Basically you would need to define the geometry first, then define analysis steps and finally invoke the analysis.\n\n"
        "Once you get used to the syntax and modelling workflow, you are more likely to write some model input files in plain text and use the '");
    suanpan_highlight("-f");
    suanpan_info("' option to directly perform the analysis, for example, '");
    suanpan_highlight("suanPan -f some_text_file");
    suanpan_info("'.\n\n"
        "To this end, a file named as 'AddAssociation.bat' (Windows) or 'suanPan.sh' (Unix) is provided to help you setup autocompletion and syntax highlighting with Sublime Text. "
        "It is placed in the same directory as the executable file. Type in '");
    suanpan_highlight("fullname");
    suanpan_info("' to see where the file is located.\n");
    restore();

    if(SUANPAN_EXIT == guide_command("fullname")) return;

    redirect();
    suanpan_info("There are other options that can be used to invoke the application. Type in '");
    suanpan_highlight("help");
    suanpan_info("' to see what options are available.\n");
    restore();

    if(SUANPAN_EXIT == guide_command("help")) return;

    redirect();
    suanpan_info("In the current interactive mode, it is also possible to load an existing model file and run it directly. "
        "Let's first echo a '");
    suanpan_highlight("benchmark");
    suanpan_info("' command to file 'benchmark.sp' using the '");
    suanpan_highlight("terminal");
    suanpan_info("' command. Type in '");
    suanpan_highlight("terminal echo benchmark >benchmark.sp");
    suanpan_info("' to do so.\n");
    restore();

    if(SUANPAN_EXIT == guide_command("terminal echo benchmark >benchmark.sp")) return;

    redirect();
    suanpan_info("Now the file can be loaded using the '");
    suanpan_highlight("file");
    suanpan_info("' command, the '");
    suanpan_highlight("benchmark");
    suanpan_info("' command we just echoed will perform some matrix solving operations, and it may take a few minutes. Type in '");
    suanpan_highlight("file benchmark.sp");
    suanpan_info("' to execute the file.\n");
    restore();

    if(SUANPAN_EXIT == guide_command("file benchmark.sp")) return;

    redirect();
    suanpan_info("In the documentation [https://tlcfem.github.io/suanPan-manual/latest/], there is an [Example] section that provides some practical examples for you to try out. "
                 "The source code repository also contains a folder named [Example] in which example usages of most models/algorithms are given. Please feel free to check that out.\n\n"
                 "Hope you will find suanPan useful. As it aims to bring the latest finite element models/algorithms to practice, you are welcome to embed your amazing research outcomes into suanPan. Type in '");
    suanpan_highlight("qrcode");
    suanpan_info("' to display a QR code for sharing. (UTF-8 is required, on Windows, some modern terminal such as Windows Terminal [https://github.com/microsoft/terminal] is recommended.)\n");
    restore();

    if(SUANPAN_EXIT == guide_command("qrcode")) return;

    redirect();
    suanpan_info("This concludes the introduction. Type '");
    suanpan_highlight("q");
    suanpan_info("' to return to the normal interactive mode.\n");
    restore();

    // ReSharper disable once CppExpressionWithoutSideEffects
    guide_command("q");
}
#pragma clang diagnostic pop

void perform_upsampling(istringstream& command) {
    string file_name;
    uword up_rate;

    if(!get_input(command, file_name, up_rate)) {
        suanpan_error("A valid file name and a valid ratio are required.\n");
        return;
    }

    string window_type = "Hamming";

    if(!get_optional_input(command, window_type)) {
        suanpan_error("A valid window type is required.\n");
        return;
    }

    mat result;

    if(is_equal(window_type, "Hamming")) result = upsampling<WindowType::Hamming>(file_name, up_rate);
    else if(is_equal(window_type, "Hann")) result = upsampling<WindowType::Hann>(file_name, up_rate);
    else if(is_equal(window_type, "Blackman")) result = upsampling<WindowType::Blackman>(file_name, up_rate);
    else if(is_equal(window_type, "BlackmanNuttall")) result = upsampling<WindowType::BlackmanNuttall>(file_name, up_rate);
    else if(is_equal(window_type, "BlackmanHarris")) result = upsampling<WindowType::BlackmanHarris>(file_name, up_rate);
    else if(is_equal(window_type, "FlatTop")) result = upsampling<WindowType::FlatTop>(file_name, up_rate);

    if(result.empty())
        suanpan_error("Fail to perform upsampling, please ensure the input is equally spaced and stored in two columns.\n");

    if(!result.save(file_name += "_upsampled", raw_ascii))
        suanpan_error("Fail to save to file.\n");
    else
        suanpan_info("Data is saved to file \"{}\".\n", file_name);
}

void perform_response_spectrum(istringstream& command) {
    string motion_name, period_name;
    if(!get_input(command, motion_name, period_name)) {
        suanpan_error("Valid file names for ground motion and period vector are required.\n");
        return;
    }

    std::error_code code;
    mat motion, period;
    if(!fs::exists(motion_name, code) || !motion.load(motion_name, raw_ascii) || motion.empty()) {
        suanpan_error("A valid ground motion stored in either one or two columns is required.\n");
        return;
    }
    if(!fs::exists(period_name, code) || !period.load(period_name, raw_ascii) || period.empty()) {
        suanpan_error("A valid period vector stored in one column is required.\n");
        return;
    }

    double interval = 0.;
    if(1llu == motion.n_cols) {
        if(!get_input(command, interval) || interval <= 0.) {
            suanpan_error("A valid sampling interval is required.\n");
            return;
        }
    }
    else {
        const vec time_diff = diff(motion.col(0));
        motion = motion.col(1);
        interval = mean(time_diff);

        if(mean(arma::abs(diff(time_diff))) > 1E-8)
            suanpan_warning("Please ensure the ground motion is equally spaced.\n");
    }

    double damping_ratio = 0.;
    if(!get_input(command, damping_ratio) || damping_ratio < 0.) {
        suanpan_error("A valid damping ratio is required.\n");
        return;
    }

    // ReSharper disable once CppTooWideScopeInitStatement
    const auto spectrum = response_spectrum<double>(damping_ratio, interval, motion, period.col(0));

    if(!spectrum.save(motion_name += "_response_spectrum", raw_ascii))
        suanpan_error("Fail to save to file.\n");
    else
        suanpan_info("Data is saved to file \"{}\".\n", motion_name);
}

void perform_sdof_response(istringstream& command) {
    string motion_name;
    if(!get_input(command, motion_name)) {
        suanpan_error("A valid file name for ground motion is required.\n");
        return;
    }

    mat motion;
    if(std::error_code code; !fs::exists(motion_name, code) || !motion.load(motion_name, raw_ascii) || motion.empty()) {
        suanpan_error("A valid ground motion stored in either one or two columns is required.\n");
        return;
    }

    double interval = 0.;
    if(1llu == motion.n_cols) {
        if(!get_input(command, interval) || interval <= 0.) {
            suanpan_error("A valid sampling interval is required.\n");
            return;
        }
    }
    else {
        const vec time_diff = diff(motion.col(0));
        motion = motion.col(1);
        interval = mean(time_diff);

        if(mean(arma::abs(diff(time_diff))) > 1E-8)
            suanpan_warning("Please make sure the ground motion is equally spaced.\n");
    }

    double freq = 0.;
    if(!get_input(command, freq) || freq <= 0.) {
        suanpan_error("A valid frequency in hertz is required.\n");
        return;
    }

    double damping_ratio = 0.;
    if(!get_input(command, damping_ratio) || damping_ratio < 0.) {
        suanpan_error("A valid damping ratio is required.\n");
        return;
    }

    // ReSharper disable once CppTooWideScopeInitStatement
    const auto response = sdof_response<double>(damping_ratio, interval, freq, motion);

    if(!response.save(motion_name += "_sdof_response", raw_ascii))
        suanpan_error("Fail to save to file.\n");
    else
        suanpan_info("Data is saved to file \"{}\".\n", motion_name);
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"
int process_command(const shared_ptr<Bead>& model, istringstream& command) {
    if(nullptr == model) return SUANPAN_SUCCESS;

    string command_id;
    if(!get_input(command, command_id)) return SUANPAN_SUCCESS;

    if(is_equal(command_id, "exit") || is_equal(command_id, "quit")) return SUANPAN_EXIT;

    if(is_equal(command_id, "overview")) {
        overview();
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "file")) {
        string file_name;
        if(!get_input(command, file_name)) {
            suanpan_error("A valid file name is required.\n");
            return SUANPAN_SUCCESS;
        }

        return process_file(model, file_name.c_str());
    }

    if(is_equal(command_id, "domain")) return create_new_domain(model, command);

    if(is_equal(command_id, "enable")) return enable_object(model, command);
    if(is_equal(command_id, "disable")) return disable_object(model, command);
    if(is_equal(command_id, "mute")) return disable_object(model, command);
    if(is_equal(command_id, "erase")) return erase_object(model, command);
    if(is_equal(command_id, "delete")) return erase_object(model, command);
    if(is_equal(command_id, "remove")) return erase_object(model, command);

    const auto& domain = get_current_domain(model);

    if(is_equal(command_id, "save")) return save_object(domain, command);
    if(is_equal(command_id, "list")) return list_object(domain, command);
    if(is_equal(command_id, "suspend")) return suspend_object(domain, command);
    if(is_equal(command_id, "protect")) return protect_object(domain, command);
    if(is_equal(command_id, "set")) return set_property(domain, command);

    if(is_equal(command_id, "amplitude")) return create_new_amplitude(domain, command);
    if(is_equal(command_id, "expression")) return create_new_expression(domain, command);
    if(is_equal(command_id, "converger")) return create_new_converger(domain, command);
    if(is_equal(command_id, "constraint")) return create_new_constraint(domain, command);
    if(is_equal(command_id, "criterion")) return create_new_criterion(domain, command);
    if(is_equal(command_id, "element")) return create_new_element(domain, command);
    if(is_equal(command_id, "hdf5recorder")) return create_new_hdf5recorder(domain, command);
    if(is_equal(command_id, "import")) return create_new_external_module(domain, command);
    if(is_equal(command_id, "initial")) return create_new_initial(domain, command);
    if(is_equal(command_id, "integrator")) return create_new_integrator(domain, command);
    if(is_equal(command_id, "load")) return create_new_load(domain, command);
    if(is_equal(command_id, "mass")) return create_new_mass(domain, command);
    if(is_equal(command_id, "material")) return create_new_material(domain, command);
    if(is_equal(command_id, "modifier")) return create_new_modifier(domain, command);
    if(is_equal(command_id, "node")) return create_new_node(domain, command);
    if(is_equal(command_id, "orientation")) return create_new_orientation(domain, command);
    if(is_equal(command_id, "plainrecorder")) return create_new_plainrecorder(domain, command);
    if(is_equal(command_id, "recorder")) return create_new_recorder(domain, command);
    if(is_equal(command_id, "section")) return create_new_section(domain, command);
    if(is_equal(command_id, "solver")) return create_new_solver(domain, command);
    if(is_equal(command_id, "step")) return create_new_step(domain, command);

    if(is_equal(command_id, "nodegroup")) return create_new_nodegroup(domain, command);
    if(is_equal(command_id, "elementgroup")) return create_new_elementgroup(domain, command);
    if(is_equal(command_id, "groupgroup")) return create_new_groupgroup(domain, command);
    if(is_equal(command_id, "generate")) return create_new_generate(domain, command);
    if(is_equal(command_id, "generatebyrule")) return create_new_generatebyrule(domain, command);
    if(is_equal(command_id, "generatebypoint")) return create_new_generatebypoint(domain, command);
    if(is_equal(command_id, "generatebyplane")) return create_new_generatebyplane(domain, command);

    auto load_handler = [&] {
        command.seekg(0);
        return create_new_load(domain, command);
    };

    if(is_equal(command_id, "acceleration")) return load_handler();
    if(is_equal(command_id, "bodyforce")) return load_handler();
    if(is_equal(command_id, "groupbodyforce")) return load_handler();
    if(is_equal(command_id, "cload")) return load_handler();
    if(is_equal(command_id, "groupcload")) return load_handler();
    if(is_equal(command_id, "lineudl2d")) return load_handler();
    if(is_equal(command_id, "lineudl3d")) return load_handler();
    if(is_equal(command_id, "disp")) return load_handler();
    if(is_equal(command_id, "displacement")) return load_handler();
    if(is_equal(command_id, "dispload")) return load_handler();
    if(is_equal(command_id, "groupdisp")) return load_handler();
    if(is_equal(command_id, "groupdisplacement")) return load_handler();
    if(is_equal(command_id, "groupdispload")) return load_handler();
    if(is_equal(command_id, "supportdisplacement")) return load_handler();
    if(is_equal(command_id, "supportvelocity")) return load_handler();
    if(is_equal(command_id, "supportacceleration")) return load_handler();

    auto constraint_handler = [&] {
        command.seekg(0);
        return create_new_constraint(domain, command);
    };

    if(is_equal(command_id, "fix")) return constraint_handler();
    if(is_equal(command_id, "penaltybc")) return constraint_handler();
    if(is_equal(command_id, "grouppenaltybc")) return constraint_handler();
    if(is_equal(command_id, "fix2")) return constraint_handler();
    if(is_equal(command_id, "multiplierbc")) return constraint_handler();
    if(is_equal(command_id, "groupmultiplierbc")) return constraint_handler();
    if(is_equal(command_id, "fixedlength2d")) return constraint_handler();
    if(is_equal(command_id, "fixedlength3d")) return constraint_handler();
    if(is_equal(command_id, "particlecollision2d")) return constraint_handler();
    if(is_equal(command_id, "particlecollision3d")) return constraint_handler();
    if(is_equal(command_id, "finiterigidwall") || is_equal(command_id, "finiterigidwallpenalty")) return constraint_handler();
    if(is_equal(command_id, "finiterigidwallmultiplier")) return constraint_handler();
    if(is_equal(command_id, "rigidwall") || is_equal(command_id, "rigidwallpenalty")) return constraint_handler();
    if(is_equal(command_id, "rigidwallmultiplier")) return constraint_handler();
    if(is_equal(command_id, "restitutionwall") || is_equal(command_id, "restitutionwallpenalty")) return constraint_handler();
    if(is_equal(command_id, "finiterestitutionwall") || is_equal(command_id, "finiterestitutionwallpenalty")) return constraint_handler();
    if(is_equal(command_id, "mpc")) return constraint_handler();

    if(is_equal(command_id, "materialtest1d")) return test_material1d(domain, command);
    if(is_equal(command_id, "materialtest2d")) return test_material2d(domain, command);
    if(is_equal(command_id, "materialtest3d")) return test_material3d(domain, command);
    if(is_equal(command_id, "materialtestwithbase3d")) return test_material_with_base3d(domain, command);
    if(is_equal(command_id, "materialtestbyload1d")) return test_material_by_load1d(domain, command);
    if(is_equal(command_id, "materialtestbyload2d")) return test_material_by_load2d(domain, command);
    if(is_equal(command_id, "materialtestbyload3d")) return test_material_by_load3d(domain, command);
    if(is_equal(command_id, "materialtestbyloadwithbase3d")) return test_material_by_load_with_base3d(domain, command);
    if(is_equal(command_id, "materialtestbystrainhistory")) return test_material_by_strain_history(domain, command);
    if(is_equal(command_id, "materialtestbystresshistory")) return test_material_by_stress_history(domain, command);

    if(is_equal(command_id, "qrcode")) {
        qrcode();
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "plot")) return vtk_parser(domain, command);

    if(is_equal(command_id, "peek")) return print_info(domain, command);

    if(is_equal(command_id, "command")) return print_command();

    if(is_equal(command_id, "example")) return run_example();

    if(is_equal(command_id, "precheck")) return model->precheck();

    if(is_equal(command_id, "analyze") || is_equal(command_id, "analyse")) {
        const auto code = model->analyze();
        suanpan_info("\n");
        return code;
    }

    if(is_equal(command_id, "fullname")) {
        suanpan_info("{}\n", SUANPAN_EXE.generic_string());
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "pwd")) {
        if(command.eof())
            suanpan_info("{}\n", fs::current_path().generic_string());
        else if(string path; get_input(command, path)) {
            std::error_code code;
            fs::current_path(path, code);
            if(0 != code.value())
                suanpan_error("Fail to set path \"{}\"\n", code.category().message(code.value()));
        }
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "help")) {
        print_helper();
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "benchmark")) return benchmark();

    if(is_equal(command_id, "clear")) {
        domain->wait();

        auto flag = true;
        for(const auto& t_integrator : domain->get_integrator_pool())
            if(t_integrator->get_domain() != nullptr) {
                t_integrator->clear_status();
                flag = false;
            }

        if(flag) domain->clear_status();

        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "reset")) {
        domain->wait();

        auto flag = true;
        for(const auto& t_integrator : domain->get_integrator_pool())
            if(t_integrator->get_domain() != nullptr) {
                t_integrator->reset_status();
                flag = false;
            }

        if(flag) domain->reset_status();

        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "summary")) {
        domain->summary();
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "terminal")) {
        execute_command(command);
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "upsampling")) {
        perform_upsampling(command);
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "response_spectrum")) {
        perform_response_spectrum(command);
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "sdof_response")) {
        perform_sdof_response(command);
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "version"))
        print_version();
    else
        suanpan_error("Command \"{}\" not found.\n", command.str());

    return SUANPAN_SUCCESS;
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"
int process_file(const shared_ptr<Bead>& model, const char* file_name) {
    std::vector<string> file_list;
    file_list.reserve(9);

    string str_name(file_name);
    file_list.emplace_back(str_name);
    file_list.emplace_back(str_name + ".supan");
    file_list.emplace_back(str_name + ".sp");

    suanpan::to_lower(str_name);
    file_list.emplace_back(str_name);
    file_list.emplace_back(str_name + ".supan");
    file_list.emplace_back(str_name + ".sp");

    suanpan::to_upper(str_name);
    file_list.emplace_back(str_name);
    file_list.emplace_back(str_name + ".SUPAN");
    file_list.emplace_back(str_name + ".SP");

    ifstream input_file;

    for(const auto& file : file_list) {
        input_file.open(fs::path(file));
        if(input_file.is_open()) break;
    }

    if(!input_file.is_open()) {
        suanpan_error("Cannot open the input file \"{}\".\n", fs::path(file_name).generic_string());
        return SUANPAN_EXIT;
    }

    string all_line, command_line;
    while(!getline(input_file, command_line).fail())
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
                if(istringstream tmp_str(all_line); process_command(model, tmp_str) == SUANPAN_EXIT) return SUANPAN_EXIT;
                all_line.clear();
            }
        }

    return SUANPAN_SUCCESS;
}
#pragma clang diagnostic pop

int create_new_domain(const shared_ptr<Bead>& model, istringstream& command) {
    unsigned domain_id;
    if(!get_input(command, domain_id)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    model->set_current_domain_tag(domain_id);

    if(auto& tmp_domain = get_domain(model, domain_id); nullptr != tmp_domain)
        suanpan_info("Switch to domain {}.\n", domain_id);
    else {
        tmp_domain = make_shared<Domain>(domain_id);
        if(nullptr != tmp_domain)
            suanpan_info("Domain {} is successfully created.\n", domain_id);
    }

    return SUANPAN_SUCCESS;
}

int disable_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(nullptr == domain) {
        suanpan_error("A valid domain is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "domain")) while(get_input(command, tag)) model->disable_domain(tag);
    else if(is_equal(object_type, "amplitude"))
        while(get_input(command, tag)) domain->disable_amplitude(tag);
    else if(is_equal(object_type, "expression"))
        while(get_input(command, tag)) domain->disable_expression(tag);
    else if(is_equal(object_type, "constraint"))
        while(get_input(command, tag)) domain->disable_constraint(tag);
    else if(is_equal(object_type, "converger"))
        while(get_input(command, tag)) domain->disable_converger(tag);
    else if(is_equal(object_type, "criterion"))
        while(get_input(command, tag)) domain->disable_criterion(tag);
    else if(is_equal(object_type, "element"))
        while(get_input(command, tag)) domain->disable_element(tag);
    else if(is_equal(object_type, "group"))
        while(get_input(command, tag)) domain->disable_group(tag);
    else if(is_equal(object_type, "integrator"))
        while(get_input(command, tag)) domain->disable_integrator(tag);
    else if(is_equal(object_type, "load"))
        while(get_input(command, tag)) domain->disable_load(tag);
    else if(is_equal(object_type, "material"))
        while(get_input(command, tag)) domain->disable_material(tag);
    else if(is_equal(object_type, "modifier"))
        while(get_input(command, tag)) domain->disable_modifier(tag);
    else if(is_equal(object_type, "node"))
        while(get_input(command, tag)) domain->disable_node(tag);
    else if(is_equal(object_type, "orientation")) while(get_input(command, tag)) domain->disable_orientation(tag);
    else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->disable_recorder(tag);
    else if(is_equal(object_type, "section")) while(get_input(command, tag)) domain->disable_section(tag);
    else if(is_equal(object_type, "solver")) while(get_input(command, tag)) domain->disable_solver(tag);
    else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->disable_step(tag);
    else if(is_equal(object_type, "print")) SUANPAN_PRINT = false;

    return SUANPAN_SUCCESS;
}

int enable_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(nullptr == domain) {
        suanpan_error("A valid domain is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "domain")) while(get_input(command, tag)) model->enable_domain(tag);
    else if(is_equal(object_type, "amplitude"))
        while(get_input(command, tag)) domain->enable_amplitude(tag);
    else if(is_equal(object_type, "expression"))
        while(get_input(command, tag)) domain->enable_expression(tag);
    else if(is_equal(object_type, "constraint"))
        while(get_input(command, tag)) domain->enable_constraint(tag);
    else if(is_equal(object_type, "converger"))
        while(get_input(command, tag)) domain->enable_converger(tag);
    else if(is_equal(object_type, "criterion"))
        while(get_input(command, tag)) domain->enable_criterion(tag);
    else if(is_equal(object_type, "element"))
        while(get_input(command, tag)) domain->enable_element(tag);
    else if(is_equal(object_type, "group"))
        while(get_input(command, tag)) domain->enable_group(tag);
    else if(is_equal(object_type, "integrator"))
        while(get_input(command, tag)) domain->enable_integrator(tag);
    else if(is_equal(object_type, "load"))
        while(get_input(command, tag)) domain->enable_load(tag);
    else if(is_equal(object_type, "material"))
        while(get_input(command, tag)) domain->enable_material(tag);
    else if(is_equal(object_type, "modifier"))
        while(get_input(command, tag)) domain->enable_modifier(tag);
    else if(is_equal(object_type, "node"))
        while(get_input(command, tag)) domain->enable_node(tag);
    else if(is_equal(object_type, "orientation")) while(get_input(command, tag)) domain->enable_orientation(tag);
    else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->enable_recorder(tag);
    else if(is_equal(object_type, "section")) while(get_input(command, tag)) domain->enable_section(tag);
    else if(is_equal(object_type, "solver")) while(get_input(command, tag)) domain->enable_solver(tag);
    else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->enable_step(tag);
    else if(is_equal(object_type, "all")) domain->enable_all();
    else if(is_equal(object_type, "print")) SUANPAN_PRINT = true;

    return SUANPAN_SUCCESS;
}

int erase_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(nullptr == domain) {
        suanpan_error("A valid domain is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "domain")) while(get_input(command, tag)) model->erase_domain(tag);
    else if(is_equal(object_type, "amplitude"))
        while(get_input(command, tag)) domain->erase_amplitude(tag);
    else if(is_equal(object_type, "expression"))
        while(get_input(command, tag)) domain->erase_expression(tag);
    else if(is_equal(object_type, "constraint"))
        while(get_input(command, tag)) domain->erase_constraint(tag);
    else if(is_equal(object_type, "converger"))
        while(get_input(command, tag)) domain->erase_converger(tag);
    else if(is_equal(object_type, "criterion"))
        while(get_input(command, tag)) domain->erase_criterion(tag);
    else if(is_equal(object_type, "element"))
        while(get_input(command, tag)) domain->erase_element(tag);
    else if(is_equal(object_type, "group"))
        while(get_input(command, tag)) domain->erase_group(tag);
    else if(is_equal(object_type, "integrator"))
        while(get_input(command, tag)) domain->erase_integrator(tag);
    else if(is_equal(object_type, "load"))
        while(get_input(command, tag)) domain->erase_load(tag);
    else if(is_equal(object_type, "material"))
        while(get_input(command, tag)) domain->erase_material(tag);
    else if(is_equal(object_type, "modifier"))
        while(get_input(command, tag)) domain->erase_modifier(tag);
    else if(is_equal(object_type, "node"))
        while(get_input(command, tag)) domain->erase_node(tag);
    else if(is_equal(object_type, "orientation")) while(get_input(command, tag)) domain->erase_orientation(tag);
    else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->erase_recorder(tag);
    else if(is_equal(object_type, "section")) while(get_input(command, tag)) domain->erase_section(tag);
    else if(is_equal(object_type, "solver")) while(get_input(command, tag)) domain->erase_solver(tag);
    else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->erase_step(tag);

    return SUANPAN_SUCCESS;
}

int save_object(const shared_ptr<DomainBase>& domain, istringstream& command) {
    if(nullptr == domain) {
        suanpan_error("A valid domain is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_id;
    if(!get_input(command, object_id)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(object_id, "Recorder")) {
        unsigned tag;
        while(get_input(command, tag)) if(domain->find_recorder(tag)) domain->get_recorder(tag)->save();
    }
    else if(is_equal(object_id, "Stiffness")) {
        string name = "K";
        if(!command.eof() && !get_input(command, name)) name = "K";
        domain->get_factory()->get_stiffness()->save(name.c_str());
    }
    else if(is_equal(object_id, "Mass")) {
        string name = "M";
        if(!command.eof() && !get_input(command, name)) name = "M";
        domain->get_factory()->get_mass()->save(name.c_str());
    }
    else if(is_equal(object_id, "Damping")) {
        string name = "C";
        if(!command.eof() && !get_input(command, name)) name = "C";
        domain->get_factory()->get_damping()->save(name.c_str());
    }
    else if(is_equal(object_id, "Model")) {
        string name = "Model.h5";
        if(!command.eof() && !get_input(command, name)) name = "Model.h5";
        domain->save(name);
    }

    return SUANPAN_SUCCESS;
}

int list_object(const shared_ptr<DomainBase>& domain, istringstream& command) {
    if(nullptr == domain) {
        suanpan_error("A valid domain is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    vector<unsigned> list;
    if(is_equal(object_type, "converger")) for(const auto& I : domain->get_converger_pool()) list.emplace_back(I->get_tag());
    else if(is_equal(object_type, "constraint")) for(const auto& I : domain->get_constraint_pool()) list.emplace_back(I->get_tag());
    else if(is_equal(object_type, "element")) for(const auto& I : domain->get_element_pool()) list.emplace_back(I->get_tag());
    else if(is_equal(object_type, "load")) for(const auto& I : domain->get_load_pool()) list.emplace_back(I->get_tag());
    else if(is_equal(object_type, "material")) for(const auto& I : domain->get_material_pool()) list.emplace_back(I->get_tag());
    else if(is_equal(object_type, "node")) for(const auto& I : domain->get_node_pool()) list.emplace_back(I->get_tag());
    else if(is_equal(object_type, "recorder")) for(const auto& I : domain->get_recorder_pool()) list.emplace_back(I->get_tag());

    suanpan_info("This domain has the following {}s:", object_type);
    for(const auto I : list)
        suanpan_info("\t{}", I);
    suanpan_info(".\n");

    return SUANPAN_SUCCESS;
}

int suspend_object(const shared_ptr<DomainBase>& domain, istringstream& command) {
    if(nullptr == domain) {
        suanpan_error("A valid domain is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    const auto step_tag = domain->get_current_step_tag();

    if(unsigned tag; is_equal(object_type, "constraint")) while(!command.eof() && get_input(command, tag)) { if(domain->find_constraint(tag)) domain->get_constraint(tag)->set_end_step(step_tag); }
    else if(is_equal(object_type, "load")) while(!command.eof() && get_input(command, tag)) { if(domain->find_load(tag)) domain->get_load(tag)->set_end_step(step_tag); }

    return SUANPAN_SUCCESS;
}

int protect_object(const shared_ptr<DomainBase>& domain, istringstream& command) {
    if(nullptr == domain) {
        suanpan_error("A valid domain is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "element")) while(!command.eof() && get_input(command, tag)) { if(domain->find<Element>(tag)) domain->get<Element>(tag)->guard(); }
    else if(is_equal(object_type, "node")) while(!command.eof() && get_input(command, tag)) { if(domain->find<Node>(tag)) domain->get<Node>(tag)->guard(); }

    return SUANPAN_SUCCESS;
}

int create_new_nodegroup(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    uword value;
    vector<uword> value_pool;
    while(get_input(command, value)) value_pool.push_back(value);

    if(!domain->insert(make_shared<NodeGroup>(tag, value_pool)))
        suanpan_error("Fail to create new node group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_elementgroup(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    uword value;
    vector<uword> value_pool;
    while(get_input(command, value)) value_pool.push_back(value);

    if(!domain->insert(make_shared<ElementGroup>(tag, value_pool)))
        suanpan_error("Fail to create new element group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_generate(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("A valid type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    int start, interval, end;
    if(!get_input(command, start)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(!get_input(command, interval)) {
        interval = 1;
        end = start;
    }
    else if(!get_input(command, end)) {
        end = interval;
        interval = end > start ? 1 : -1;
    }

    if(0 == interval) interval = 1;

    if(start == end) interval = 1;
    else if(start < end && interval < 0 || start > end && interval > 0) interval = -interval;

    vector<uword> tag_pool;

    tag_pool.reserve(std::max(1, (end - start) / interval + 1));

    while(start <= end) {
        tag_pool.emplace_back(start);
        start += interval;
    }

    if(is_equal(type, "nodegroup") && !domain->insert(make_shared<NodeGroup>(tag, tag_pool)))
        suanpan_error("Fail to create new node group.\n");
    else if(is_equal(type, "elementgroup") && !domain->insert(make_shared<ElementGroup>(tag, tag_pool)))
        suanpan_error("Fail to create new element group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_generatebyrule(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("A valid type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof;
    if(!get_input(command, dof)) {
        suanpan_error("A valid dof identifier is required.\n");
        return SUANPAN_SUCCESS;
    }

    double para;
    vector<double> pool;
    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(is_equal(type, "nodegroup") && !domain->insert(make_shared<NodeGroup>(tag, dof, pool)))
        suanpan_error("Fail to create new node group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_generatebyplane(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("A valid type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    double para;
    vector<double> pool;
    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(pool.empty()) return SUANPAN_SUCCESS;

    if(is_equal(type, "nodegroup") && !domain->insert(make_shared<NodeGroup>(tag, pool)))
        suanpan_error("Fail to create new node group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_generatebypoint(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("A valid type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    double para;
    vector<double> pool;
    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(pool.size() % 2 == 0) {
        if(const auto size = static_cast<long long>(pool.size()) / 2; is_equal(type, "nodegroup") && !domain->insert(make_shared<NodeGroup>(tag, vector(pool.begin(), pool.begin() + size), vector(pool.end() - size, pool.end()))))
            suanpan_error("Fail to create new node group.\n");
    }

    return SUANPAN_SUCCESS;
}

int create_new_groupgroup(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    uword para;
    vector<uword> pool;
    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(!domain->insert(make_shared<GroupGroup>(tag, pool)))
        suanpan_error("Fail to create new group of groups.\n");

    return SUANPAN_SUCCESS;
}

int create_new_external_module(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string library_name;

    if(!get_input(command, library_name)) {
        suanpan_error("A valid module name is required.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    for(const auto& I : domain->get_external_module_pool())
        if(is_equal(I->library_name, library_name) || I->locate_cpp_module(library_name)) {
            code = 1;
            break;
        }

    if(0 == code) domain->insert(make_shared<ExternalModule>(library_name));

    return SUANPAN_SUCCESS;
}

int create_new_initial(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string variable_type;
    if(!get_input(command, variable_type)) {
        suanpan_error("A valid variable type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal("material", variable_type)) {
        string state_type;
        if(!get_input(command, state_type)) {
            suanpan_error("A valid state type is required.\n");
            return SUANPAN_SUCCESS;
        }

        unsigned mat_tag;
        if(!get_input(command, mat_tag)) {
            suanpan_error("A valid material tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        vector<double> para;
        while(!command.eof()) if(double input; get_input(command, input)) para.emplace_back(input);

        if(is_equal("history", state_type) && domain->find_material(mat_tag)) domain->get_material(mat_tag)->set_initial_history(para);

        return SUANPAN_SUCCESS;
    }
    if(is_equal("angularvelocity", variable_type) || is_equal("avel", variable_type)) {
        vec magnitude(3);
        for(auto& I : magnitude)
            if(!get_input(command, I)) {
                suanpan_error("A valid magnitude is required.\n");
                return SUANPAN_SUCCESS;
            }

        unsigned ref_node;
        if(!get_input(command, ref_node) || !domain->find_node(ref_node)) {
            suanpan_error("A valid reference node tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        auto& t_ref_node = domain->get_node(ref_node);
        auto t_ref_coor = t_ref_node->get_coordinate();
        t_ref_coor.resize(3);

        while(!command.eof()) {
            if(unsigned node_tag; get_input(command, node_tag) && domain->find_node(node_tag)) {
                auto& t_node = domain->get_node(node_tag);
                auto t_coor = t_node->get_coordinate();
                t_coor.resize(3);
                t_node->update_current_velocity(cross(magnitude, t_coor - t_ref_coor));
            }
            else {
                suanpan_error("A valid node tag is required.\n");
                return SUANPAN_SUCCESS;
            }
        }

        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("A valid magnitude is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_tag;
    if(!get_input(command, dof_tag)) {
        suanpan_error("A valid dof identifier is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal("displacement", variable_type) || is_equal("disp", variable_type))
        while(!command.eof())
            if(unsigned node_tag; get_input(command, node_tag) && domain->find_node(node_tag)) {
                auto& t_node = domain->get_node(node_tag);
                auto t_variable = t_node->get_current_displacement();
                if(t_variable.n_elem < dof_tag) t_variable.resize(dof_tag);
                t_variable(dof_tag - 1) = magnitude;
                t_node->update_current_displacement(t_variable);
            }
            else {
                suanpan_error("A valid node tag is required.\n");
                return SUANPAN_SUCCESS;
            }
    else if(is_equal("velocity", variable_type) || is_equal("vel", variable_type))
        while(!command.eof())
            if(unsigned node_tag; get_input(command, node_tag) && domain->find_node(node_tag)) {
                auto& t_node = domain->get_node(node_tag);
                auto t_variable = t_node->get_current_velocity();
                if(t_variable.n_elem < dof_tag) t_variable.resize(dof_tag);
                t_variable(dof_tag - 1) = magnitude;
                t_node->update_current_velocity(t_variable);
            }
            else {
                suanpan_error("A valid node tag is required.\n");
                return SUANPAN_SUCCESS;
            }
    else if(is_equal("acceleration", variable_type) || is_equal("acc", variable_type))
        while(!command.eof()) {
            if(unsigned node_tag; get_input(command, node_tag) && domain->find_node(node_tag)) {
                auto& t_node = domain->get_node(node_tag);
                auto t_variable = t_node->get_current_acceleration();
                if(t_variable.n_elem < dof_tag) t_variable.resize(dof_tag);
                t_variable(dof_tag - 1) = magnitude;
                t_node->update_current_acceleration(t_variable);
            }
            else {
                suanpan_error("A valid node tag is required.\n");
                return SUANPAN_SUCCESS;
            }
        }

    return SUANPAN_SUCCESS;
}

int create_new_node(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned node_id;
    if(!get_input(command, node_id)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vector<double> coor;
    double X;
    while(get_input(command, X)) coor.push_back(X);

    if(!domain->insert(make_shared<Node>(node_id, vec(coor))))
        suanpan_error("Fail to create new node via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int set_property(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string property_id;
    if(!get_input(command, property_id)) {
        suanpan_error("A valid property type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(property_id, "output_folder")) {
        string value;

        if(!get_input(command, value)) {
            suanpan_error("A valid value is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(value, "$pwd")) SUANPAN_OUTPUT = canonical(fs::current_path());
        else {
            fs::path new_path = value;
            if(new_path.is_relative()) new_path = SUANPAN_OUTPUT / new_path;

            if(!exists(new_path)) {
                std::error_code code;
                create_directories(new_path, code);
                if(0 != code.value()) {
                    suanpan_error("Cannot create \"{}\".\n", new_path.generic_string());
                    return SUANPAN_SUCCESS;
                }
            }

            SUANPAN_OUTPUT = canonical(new_path);
        }

        suanpan_info("{}\n", SUANPAN_OUTPUT.generic_string());
        return SUANPAN_SUCCESS;
    }
    if(is_equal(property_id, "num_threads")) {
        if(int value; get_input(command, value)) SUANPAN_NUM_THREADS = value;
        else
            suanpan_error("A valid value is required.\n");

        return SUANPAN_SUCCESS;
    }
    if(is_equal(property_id, "screen_output")) {
        if(string value; get_input(command, value)) SUANPAN_PRINT = is_true(value);
        else
            suanpan_error("A valid value is required.\n");

        return SUANPAN_SUCCESS;
    }

    if(is_equal(property_id, "color_model")) {
        if(string value; !get_input(command, value))
            suanpan_error("A valid value is required.\n");
        else if(is_equal("WP", value)) domain->set_color_model(ColorMethod::WP);
        else if(is_equal("MIS", value)) domain->set_color_model(ColorMethod::MIS);
        else domain->set_color_model(ColorMethod::OFF);

        return SUANPAN_SUCCESS;
    }
    if(is_equal(property_id, "constraint_multiplier")) {
        double value;
        get_input(command, value) ? set_constraint_multiplier(value) : suanpan_error("A valid value is required.\n");

        return SUANPAN_SUCCESS;
    }
    if(is_equal(property_id, "load_multiplier")) {
        double value;
        get_input(command, value) ? set_load_multiplier(value) : suanpan_error("A valid value is required.\n");

        return SUANPAN_SUCCESS;
    }

    if(is_equal(property_id, "linear_system")) {
        domain->set_attribute(ModalAttribute::LinearSystem);

        return SUANPAN_SUCCESS;
    }

    if(domain->get_current_step_tag() == 0) return SUANPAN_SUCCESS;

    const auto& t_step = domain->get_current_step();

    if(is_equal(property_id, "fixed_step_size")) {
        string value;
        get_input(command, value) ? t_step->set_fixed_step_size(is_true(value)) : suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "symm_mat")) {
        string value;
        get_input(command, value) ? t_step->set_symm(is_true(value)) : suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "band_mat")) {
        string value;
        get_input(command, value) ? t_step->set_band(is_true(value)) : suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "sparse_mat")) {
        string value;
        get_input(command, value) ? t_step->set_sparse(is_true(value)) : suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "iterative_refinement")) {
        if(unsigned value; get_input(command, value)) t_step->set_refinement(value);
        else
            suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "system_solver")) {
        if(string value; !get_input(command, value))
            suanpan_error("A valid value is required.\n");
        else if(is_equal(value, "LAPACK")) t_step->set_system_solver(SolverType::LAPACK);
        else if(is_equal(value, "SPIKE")) t_step->set_system_solver(SolverType::SPIKE);
        else if(is_equal(value, "SUPERLU")) t_step->set_system_solver(SolverType::SUPERLU);
        else if(is_equal(value, "MUMPS")) t_step->set_system_solver(SolverType::MUMPS);
#ifdef SUANPAN_CUDA
        else if(is_equal(value, "CUDA")) t_step->set_system_solver(SolverType::CUDA);
#endif
#ifdef SUANPAN_MKL
        else if(is_equal(value, "PARDISO")) t_step->set_system_solver(SolverType::PARDISO);
        else if(is_equal(value, "FGMRES")) t_step->set_system_solver(SolverType::FGMRES);
#endif
        else if(is_equal(value, "GMRES")) t_step->set_system_solver(IterativeSolver::GMRES);
        else if(is_equal(value, "BICGSTAB")) t_step->set_system_solver(IterativeSolver::BICGSTAB);
        else if(is_equal(value, "NONE")) t_step->set_system_solver(IterativeSolver::NONE);
        else
            suanpan_error("A valid solver type is required.\n");
    }
    else if(is_equal(property_id, "preconditioner")) {
        if(string value; !get_input(command, value))
            suanpan_error("A valid value is required.\n");
        else if(is_equal(value, "NONE")) t_step->set_preconditioner(PreconditionerType::NONE);
        else if(is_equal(value, "JACOBI")) t_step->set_preconditioner(PreconditionerType::JACOBI);
#ifndef SUANPAN_SUPERLUMT
        else if(is_equal(value, "ILU")) t_step->set_preconditioner(PreconditionerType::ILU);
#endif
        else
            suanpan_error("A valid solver type is required.\n");
    }
    else if(is_equal(property_id, "precision")) {
        if(string value; !get_input(command, value))
            suanpan_error("A valid value is required.\n");
        else if(is_equal(value, "DOUBLE") || is_equal(value, "FULL")) t_step->set_precision(Precision::FULL);
        else if(is_equal(value, "SINGLE") || is_equal(value, "MIXED")) t_step->set_precision(Precision::MIXED);
        else
            suanpan_error("A valid precision is required.\n");
    }
    else if(is_equal(property_id, "tolerance")) {
        double value;
        get_input(command, value) ? t_step->set_tolerance(value) : suanpan_error("A valid value is required.\n");
    }
#ifdef SUANPAN_MKL
    else if(is_equal(property_id, "fgmres_tolerance")) {
        double value;
        get_input(command, value) ? t_step->set_tolerance(value) : suanpan_error("A valid value is required.\n");
    }
#endif
    else if(is_equal(property_id, "ini_step_size")) {
        double step_time;
        get_input(command, step_time) ? t_step->set_ini_step_size(step_time) : suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "min_step_size")) {
        double step_time;
        get_input(command, step_time) ? t_step->set_min_step_size(step_time) : suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "max_step_size")) {
        double step_time;
        get_input(command, step_time) ? t_step->set_max_step_size(step_time) : suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "max_iteration")) {
        unsigned max_number;
        get_input(command, max_number) ? t_step->set_max_substep(max_number) : suanpan_error("A valid value is required.\n");
    }
    else if(is_equal(property_id, "eigen_number")) {
        if(unsigned eigen_number; get_input(command, eigen_number)) {
            if(const auto eigen_step = std::dynamic_pointer_cast<Frequency>(t_step); nullptr == eigen_step)
                suanpan_error("Cannot set eigen number for non-eigen step.\n");
            else eigen_step->set_eigen_number(eigen_number);
        }
        else
            suanpan_error("A valid eigen number is required.\n");
    }

    return SUANPAN_SUCCESS;
}

int print_info(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "node"))
        while(get_input(command, tag)) {
            if(domain->find_node(tag)) {
                get_node(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "element"))
        while(get_input(command, tag)) {
            if(domain->find_element(tag)) {
                get_element(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "material"))
        while(get_input(command, tag)) {
            if(domain->find_material(tag)) {
                get_material(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "constraint"))
        while(get_input(command, tag)) {
            if(domain->find_constraint(tag)) {
                get_constraint(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "recorder"))
        while(get_input(command, tag)) {
            if(domain->find_recorder(tag)) {
                get_recorder(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "solver"))
        while(get_input(command, tag)) {
            if(domain->find_solver(tag)) {
                get_solver(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "integrator"))
        while(get_input(command, tag)) {
            if(domain->find_integrator(tag)) {
                get_integrator(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "group"))
        while(get_input(command, tag)) {
            if(domain->find_group(tag)) {
                get_group(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "nodegroup"))
        while(get_input(command, tag)) {
            if(domain->find_group(tag))
                for(const auto t_node : get_group(domain, tag)->get_pool())
                    if(domain->find<Node>(t_node)) {
                        get_node(domain, static_cast<unsigned>(t_node))->print();
                        suanpan_info("\n");
                    }
        }
    else if(is_equal(object_type, "amplitude"))
        while(get_input(command, tag)) {
            if(domain->find_amplitude(tag)) {
                get_amplitude(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "expression"))
        while(get_input(command, tag)) {
            if(domain->find_expression(tag)) {
                get_expression(domain, tag)->print();
                suanpan_info("\n");
            }
        }
    else if(is_equal(object_type, "eigenvalue")) {
        domain->get_factory()->get_eigenvalue().print("Eigenvalues:");
        suanpan_info("\n");
    }
    else if(is_equal(object_type, "output_folder"))
        suanpan_info("{}\n", SUANPAN_OUTPUT.generic_string());

    return SUANPAN_SUCCESS;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"
int run_example() {
    const auto new_model = make_shared<Bead>();

    suanpan_info("====================================================\n");
    suanpan_info("-> A Minimum Example: Elastic Truss Under Tension <-\n");
    suanpan_info("====================================================\n");

    constexpr auto wait_time = 1000;

    auto run_command = [&](const string& command_string) {
        suanpan_highlight("\t{}\n", command_string);
        auto command = istringstream(command_string);
        process_command(new_model, command);
        std::this_thread::sleep_for(std::chrono::milliseconds(wait_time));
    };

    // node
    suanpan_info("--> create two nodes at (0,0) and (2,0):\n");
    run_command("node 1 0 0");
    run_command("node 2 2 0");

    // material
    suanpan_info("--> create material model (elastic modulus 52):\n");
    run_command("material Elastic1D 1 52");

    // element
    suanpan_info("--> create a truss element connecting nodes 1 and 2:\n");
    run_command("element T2D2 1 1 2 1 93");

    // boundary condition and load
    suanpan_info("--> define boundary condition and load:\n");
    run_command("fix 1 1 1");
    run_command("fix 2 2 1 2");
    run_command("displacement 1 0 1.4 1 2");

    // step
    suanpan_info("--> define a static step:\n");
    run_command("step static 1");

    // analyze
    suanpan_info("--> perform the analysis:\n");
    run_command("analyze");

    // analyze
    suanpan_info("--> check nodal force (P=UEA/L=1.4*52*93/2=3385.2):\n");
    run_command("peek node 2");

    // clean up
    suanpan_info("--> clean up and it's your turn!\n");

    suanpan_info("====================================================\n");
    return SUANPAN_SUCCESS;
}
#pragma clang diagnostic pop

int print_command() {
    suanpan_info("The available commands are listed. Please check online manual for reference. https://tlcfem.gitbook.io/suanpan-manual/\n");

    constexpr auto format = "    {:<30}  {}\n";
    suanpan_info(format, "acceleration", "define acceleration");
    suanpan_info(format, "amplitude", "define amplitude");
    suanpan_info(format, "analyze/analyse", "analyse the model");
    suanpan_info(format, "benchmark", "benchmark the platform for comparison");
    suanpan_info(format, "clear", "clear model");
    suanpan_info(format, "cload", "define concentrated load");
    suanpan_info(format, "command", "list all commands");
    suanpan_info(format, "converger", "define converger");
    suanpan_info(format, "criterion", "define stopping criterion");
    suanpan_info(format, "delete/erase/remove", "delete objects");
    suanpan_info(format, "disable/mute", "disable objects");
    suanpan_info(format, "disp/displacement/dispload", "define displacement load");
    suanpan_info(format, "domain", "create/switch to domains");
    suanpan_info(format, "element", "define element");
    suanpan_info(format, "elementgroup", "define group containing element tags");
    suanpan_info(format, "enable", "enable objects");
    suanpan_info(format, "example", "establish adn execute a minimum example");
    suanpan_info(format, "exit/quit", "exit program");
    suanpan_info(format, "file", "load external files");
    suanpan_info(format, "finiterigidwall", "define rigid wall constraint with finite dimensions");
    suanpan_info(format, "fix/penaltybc", "define boundary conditions by penalty method");
    suanpan_info(format, "fix2/multiplierbc", "define boundary conditions by multiplier method");
    suanpan_info(format, "fullname", "print the full path of the program");
    suanpan_info(format, "generate", "generate node or element group by fixed interval");
    suanpan_info(format, "generatebyplane", "generate node or element group by plane");
    suanpan_info(format, "generatebypoint", "generate node or element group by line segment");
    suanpan_info(format, "generatebyrule", "generate node or element group by polynomial");
    suanpan_info(format, "groupcload", "define concentrated load based on given group");
    suanpan_info(format, "groupdisp", "define displacement load based on given group");
    suanpan_info(format, "groupmultiplierbc", "define boundary conditions by multiplier method based on given group");
    suanpan_info(format, "grouppenaltybc", "define boundary conditions by penalty method based on given group");
    suanpan_info(format, "hdf5recorder", "define recorder using hdf5 format");
    suanpan_info(format, "help", "print help information");
    suanpan_info(format, "import", "import external module");
    suanpan_info(format, "initial", "define initial condition");
    suanpan_info(format, "integrator", "define time integration algorithm");
    suanpan_info(format, "list", "list objects in the current domain");
    suanpan_info(format, "mass", "define point mass");
    suanpan_info(format, "material", "define material");
    suanpan_info(format, "materialtest1d", "test independent material modal using displacement input");
    suanpan_info(format, "materialtest2d", "test independent material modal using displacement input");
    suanpan_info(format, "materialtest3d", "test independent material modal using displacement input");
    suanpan_info(format, "materialtestbyload1d", "test independent material modal using force input");
    suanpan_info(format, "materialtestbyload2d", "test independent material modal using force input");
    suanpan_info(format, "materialtestbyload3d", "test independent material modal using force input");
    suanpan_info(format, "materialtestbyloadwithbase3d", "test independent material modal using force input");
    suanpan_info(format, "materialtestwithbase3d", "test independent material modal using displacement input");
    suanpan_info(format, "modifier", "define modifier that modifies existing modal properties");
    suanpan_info(format, "mpc", "define multi-point constraint");
    suanpan_info(format, "node", "define node");
    suanpan_info(format, "nodegroup", "define group containing node tags");
    suanpan_info(format, "orientation", "define beam section orientation");
    suanpan_info(format, "particlecollision", "define collision constraint between particles");
    suanpan_info(format, "peek", "peek current information of target object");
    suanpan_info(format, "plainrecorder", "define recorder using plain text format");
    suanpan_info(format, "plot", "plot and optionally save model");
    suanpan_info(format, "precheck", "check the model without analyse");
    suanpan_info(format, "protect", "protect objects from being disabled");
    suanpan_info(format, "pwd", "print/change current working folder");
    suanpan_info(format, "qrcode", "print a qr code");
    suanpan_info(format, "recorder", "define recorder");
    suanpan_info(format, "reset", "reset model to the previously converged state");
    suanpan_info(format, "rigidwall", "define rigid wall constraint with infinite dimensions");
    suanpan_info(format, "save", "save objects");
    suanpan_info(format, "section", "define section");
    suanpan_info(format, "set", "set properties of analysis");
    suanpan_info(format, "solver", "define solver");
    suanpan_info(format, "step", "define step");
    suanpan_info(format, "summary", "print summary for the current domain");
    suanpan_info(format, "suspend", "suspend object in current step");
    suanpan_info(format, "terminal", "execute command in terminal");
    suanpan_info(format, "version", "print version information");

    return SUANPAN_SUCCESS;
}

int execute_command(istringstream& command) {
#ifdef SUANPAN_MSVC
    std::wstringstream terminal_command;
    terminal_command << command.str().substr(command.tellg()).c_str();
    const auto code = _wsystem(terminal_command.str().c_str());
#else
    std::stringstream terminal_command;
    terminal_command << command.str().substr(command.tellg()).c_str();
    const auto code = system(terminal_command.str().c_str());
#endif

    return code;
}
