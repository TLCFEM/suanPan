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

// ReSharper disable StringLiteralTypo
// ReSharper disable IdentifierTypo
#include "command.h"

#include <Constraint/Constraint.h>
#include <Constraint/ConstraintParser.h>
#include <Constraint/Criterion/Criterion.h>
#include <Converger/Converger.h>
#include <Converger/ConvergerParser.h>
#include <Domain/Domain.h>
#include <Domain/ExternalModule.h>
#include <Domain/Group/Group.h>
#include <Domain/Group/GroupParser.h>
#include <Domain/Node.h>
#include <Element/Element.h>
#include <Element/ElementParser.h>
#include <Element/Modifier/Modifier.h>
#include <Element/Utility/Orientation.h>
#include <Element/Visualisation/vtkParser.h>
#include <Load/Amplitude/Amplitude.h>
#include <Load/Load.h>
#include <Load/LoadParser.h>
#include <Material/Material.h>
#include <Material/MaterialParser.h>
#include <Material/MaterialTester.h>
#include <Recorder/Recorder.h>
#include <Recorder/RecorderParser.h>
#include <Section/Section.h>
#include <Section/SectionParser.h>
#include <Section/SectionTester.h>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Solver.h>
#include <Solver/SolverParser.h>
#include <Step/Bead.h>
#include <Step/Frequency.h>
#include <Step/Step.h>
#include <Step/StepParser.h>
#include <Toolbox/Expression.h>
#include <Toolbox/ExpressionParser.h>
#include <Toolbox/argument.h>
#include <Toolbox/resampling.h>
#include <Toolbox/response_spectrum.h>
#include <Toolbox/thread_pool.hpp>
#include <thread>

using std::ifstream;
using std::ofstream;
using std::vector;

unsigned SUANPAN_WARNING_COUNT = 0;
unsigned SUANPAN_ERROR_COUNT = 0;
int SUANPAN_NUM_THREADS = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
int SUANPAN_NUM_NODES = comm_size;
fs::path SUANPAN_OUTPUT = fs::current_path();
extern fs::path SUANPAN_EXE;

namespace {
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

    int qrcode() {
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

        return SUANPAN_SUCCESS;
    }

    void overview() {
        const auto new_model = std::make_shared<Bead>();

        auto guide_command = [&](const std::string& target) {
            while(true) {
                std::string command_line;
                suanpan_highlight("overview -> ");
                getline(std::cin, command_line);
                if(is_equal(command_line, "q")) {
                    suanpan_info("Bye. Happy modelling!\n\n");
                    return SUANPAN_EXIT;
                }
                if(is_equal(command_line, target)) {
                    std::istringstream tmp_str(command_line);
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
            buffer_tmp = std::make_unique<std::stringbuf>();
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
        suanpan_info("Welcome to suanPan, a parallel finite element analysis software.\n\nIt can be used to solve various types of solid mechanics problems. It also provides a flexible platform for researchers to develop new numerical models and algorithms in a hassle-free manner.\n\nThe following is a quick introduction of the application, regarding its functionalities and basic usage. At any time, type in '");
        suanpan_highlight("q");
        suanpan_info("' to quit this introduction. First of all, type in '");
        suanpan_highlight("version");
        suanpan_info("' to check the license and version information.\n");
        restore();

        if(SUANPAN_EXIT == guide_command("version")) return;

        redirect();
        suanpan_info("By invoking the application without input file, you are currently in the interactive mode, which allows you to create finite element models interactively. A numerical model typically consists of nodes, elements connecting nodes, materials attached to elements, loads applied to nodes, boundary conditions imposed on nodes, etc. To analyze, a number of analysis steps need to be defined with the proper convergence criteria and solvers.\n\nType in '");
        suanpan_highlight("example");
        suanpan_info("' to check out a simple elastic cantilever example.\n");
        restore();

        if(SUANPAN_EXIT == guide_command("example")) return;

        redirect();
        suanpan_info("If you have some experience with ABAQUS scripting, you may have realised that the structure resembles that of ABAQUS *.inp files. Basically you would need to define the geometry first, then define analysis steps and finally invoke the analysis.\n\nOnce you get used to the syntax and modelling workflow, you are more likely to write some model input files in plain text and use the '");
        suanpan_highlight("-f");
        suanpan_info("' option to directly perform the analysis, for example, '");
        suanpan_highlight("suanPan -f some_text_file");
        suanpan_info("'.\n\nTo this end, a file named as 'AddAssociation.bat' (Windows) or 'suanPan.sh' (Unix) is provided to help you setup autocompletion and syntax highlighting with Sublime Text. It is placed in the same directory as the executable file. Type in '");
        suanpan_highlight("fullname");
        suanpan_info("' to see where the file is located.\n");
        restore();

        if(SUANPAN_EXIT == guide_command("fullname")) return;

        redirect();
        suanpan_info("There are a number of top-level commands. Type in '");
        suanpan_highlight("command");
        suanpan_info("' to see a list.\n");
        restore();

        if(SUANPAN_EXIT == guide_command("command")) return;

        redirect();
        suanpan_info("In the current interactive mode, it is also possible to load an existing model file and run it directly. Let's first echo a '");
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
        suanpan_info("In the documentation [https://tlcfem.github.io/suanPan-manual/latest/], there is an [Example] section that provides some practical examples for you to try out. The source code repository also contains a folder named [Example] in which example usages of most models/algorithms are given. Please feel free to check that out.\n\nHope you will find suanPan useful. As it aims to bring the latest finite element models/algorithms to practice, you are welcome to embed your amazing research outcomes into suanPan. Type in '");
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

    void perform_upsampling(std::istringstream& command) {
        std::string file_name;
        uword up_rate;

        if(!get_input(command, file_name, up_rate)) {
            suanpan_error("A valid file name and a valid ratio are required.\n");
            return;
        }

        std::string window_type = "Hamming";
        if(!get_optional_input(command, window_type)) {
            suanpan_error("A valid window type is required.\n");
            return;
        }

        auto window_size = 8llu;
        if(!get_optional_input(command, window_size)) {
            suanpan_error("A valid window size is required.\n");
            return;
        }

        const mat result = upsampling(window_type, file_name, up_rate, window_size);

        if(result.empty())
            suanpan_error("Fail to perform upsampling, please ensure the input is equally spaced and stored in two columns.\n");

        if(!result.save(file_name += "_upsampled", raw_ascii))
            suanpan_error("Fail to save to file.\n");
        else
            suanpan_info("Data is saved to file \"{}\".\n", file_name);
    }

    void perform_response_spectrum(std::istringstream& command) {
        std::string motion_name, period_name;
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

        auto interval = 0.;
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

        auto damping_ratio = 0.;
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

    void perform_sdof_response(std::istringstream& command) {
        std::string motion_name;
        if(!get_input(command, motion_name)) {
            suanpan_error("A valid file name for ground motion is required.\n");
            return;
        }

        mat motion;
        if(std::error_code code; !fs::exists(motion_name, code) || !motion.load(motion_name, raw_ascii) || motion.empty()) {
            suanpan_error("A valid ground motion stored in either one or two columns is required.\n");
            return;
        }

        auto interval = 0.;
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

        auto freq = 0.;
        if(!get_input(command, freq) || freq <= 0.) {
            suanpan_error("A valid frequency in hertz is required.\n");
            return;
        }

        auto damping_ratio = 0.;
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

    int create_new_domain(const shared_ptr<Bead>& model, std::istringstream& command) {
        unsigned domain_id;
        if(!get_input(command, domain_id)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        model->set_current_domain_tag(domain_id);

        if(auto& tmp_domain = get_domain(model, domain_id); nullptr != tmp_domain)
            suanpan_info("Switch to domain {}.\n", domain_id);
        else {
            tmp_domain = std::make_shared<Domain>(domain_id);
            if(nullptr != tmp_domain)
                suanpan_info("Domain {} is successfully created.\n", domain_id);
        }

        return SUANPAN_SUCCESS;
    }

    int create_new_external_module(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        std::string library_name;

        if(!get_input(command, library_name)) {
            suanpan_error("A valid module name is required.\n");
            return SUANPAN_SUCCESS;
        }

        auto code = 0;
        for(auto& I : domain->get_external_module_pool())
            if(is_equal(I->library_name, library_name) || I->locate_cpp_module(library_name)) {
                code = 1;
                break;
            }

        if(0 == code) domain->insert(std::make_shared<ExternalModule>(library_name));

        return SUANPAN_SUCCESS;
    }

    int create_new_initial(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        std::string variable_type;
        if(!get_input(command, variable_type)) {
            suanpan_error("A valid variable type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal("material", variable_type)) {
            std::string state_type;
            if(!get_input(command, state_type)) {
                suanpan_error("A valid state type is required.\n");
                return SUANPAN_SUCCESS;
            }

            unsigned mat_tag;
            if(!get_input(command, mat_tag)) {
                suanpan_error("A valid material tag is required.\n");
                return SUANPAN_SUCCESS;
            }

            std::vector<double> para;
            while(!command.eof())
                if(double input; get_input(command, input)) para.emplace_back(input);

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

    int create_new_node(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        unsigned node_id;
        if(!get_input(command, node_id)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::vector<double> coor;
        double X;
        while(get_input(command, X)) coor.push_back(X);

        if(!domain->insert(std::make_shared<Node>(node_id, vec(coor))))
            suanpan_error("Fail to create new node via \"{}\".\n", command.str());

        return SUANPAN_SUCCESS;
    }

    int disable_object(const shared_ptr<Bead>& model, std::istringstream& command) {
        const auto& domain = get_current_domain(model);
        if(nullptr == domain) {
            suanpan_error("A valid domain is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string object_type;
        if(!get_input(command, object_type)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(unsigned tag; is_equal(object_type, "amplitude"))
            while(get_input(command, tag)) domain->disable_amplitude(tag);
        else if(is_equal(object_type, "constraint"))
            while(get_input(command, tag)) domain->disable_constraint(tag);
        else if(is_equal(object_type, "converger"))
            while(get_input(command, tag)) domain->disable_converger(tag);
        else if(is_equal(object_type, "criterion"))
            while(get_input(command, tag)) domain->disable_criterion(tag);
        else if(is_equal(object_type, "domain"))
            while(get_input(command, tag)) model->disable_domain(tag);
        else if(is_equal(object_type, "element"))
            while(get_input(command, tag)) domain->disable_element(tag);
        else if(is_equal(object_type, "expression"))
            while(get_input(command, tag)) domain->disable_expression(tag);
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
        else if(is_equal(object_type, "orientation"))
            while(get_input(command, tag)) domain->disable_orientation(tag);
        else if(is_equal(object_type, "recorder"))
            while(get_input(command, tag)) domain->disable_recorder(tag);
        else if(is_equal(object_type, "section"))
            while(get_input(command, tag)) domain->disable_section(tag);
        else if(is_equal(object_type, "solver"))
            while(get_input(command, tag)) domain->disable_solver(tag);
        else if(is_equal(object_type, "step"))
            while(get_input(command, tag)) domain->disable_step(tag);
        else if(is_equal(object_type, "print")) SUANPAN_PRINT = false;

        return SUANPAN_SUCCESS;
    }

    int enable_object(const shared_ptr<Bead>& model, std::istringstream& command) {
        const auto& domain = get_current_domain(model);
        if(nullptr == domain) {
            suanpan_error("A valid domain is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string object_type;
        if(!get_input(command, object_type)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(unsigned tag; is_equal(object_type, "amplitude"))
            while(get_input(command, tag)) domain->enable_amplitude(tag);
        else if(is_equal(object_type, "constraint"))
            while(get_input(command, tag)) domain->enable_constraint(tag);
        else if(is_equal(object_type, "converger"))
            while(get_input(command, tag)) domain->enable_converger(tag);
        else if(is_equal(object_type, "criterion"))
            while(get_input(command, tag)) domain->enable_criterion(tag);
        else if(is_equal(object_type, "domain"))
            while(get_input(command, tag)) model->enable_domain(tag);
        else if(is_equal(object_type, "element"))
            while(get_input(command, tag)) domain->enable_element(tag);
        else if(is_equal(object_type, "expression"))
            while(get_input(command, tag)) domain->enable_expression(tag);
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
        else if(is_equal(object_type, "orientation"))
            while(get_input(command, tag)) domain->enable_orientation(tag);
        else if(is_equal(object_type, "recorder"))
            while(get_input(command, tag)) domain->enable_recorder(tag);
        else if(is_equal(object_type, "section"))
            while(get_input(command, tag)) domain->enable_section(tag);
        else if(is_equal(object_type, "solver"))
            while(get_input(command, tag)) domain->enable_solver(tag);
        else if(is_equal(object_type, "step"))
            while(get_input(command, tag)) domain->enable_step(tag);
        else if(is_equal(object_type, "all")) domain->enable_all();
        else if(is_equal(object_type, "print")) SUANPAN_PRINT = true;

        return SUANPAN_SUCCESS;
    }

    int erase_object(const shared_ptr<Bead>& model, std::istringstream& command) {
        const auto& domain = get_current_domain(model);
        if(nullptr == domain) {
            suanpan_error("A valid domain is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string object_type;
        if(!get_input(command, object_type)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(unsigned tag; is_equal(object_type, "amplitude"))
            while(get_input(command, tag)) domain->erase_amplitude(tag);
        else if(is_equal(object_type, "constraint"))
            while(get_input(command, tag)) domain->erase_constraint(tag);
        else if(is_equal(object_type, "converger"))
            while(get_input(command, tag)) domain->erase_converger(tag);
        else if(is_equal(object_type, "criterion"))
            while(get_input(command, tag)) domain->erase_criterion(tag);
        else if(is_equal(object_type, "domain"))
            while(get_input(command, tag)) model->erase_domain(tag);
        else if(is_equal(object_type, "element"))
            while(get_input(command, tag)) domain->erase_element(tag);
        else if(is_equal(object_type, "expression"))
            while(get_input(command, tag)) domain->erase_expression(tag);
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
        else if(is_equal(object_type, "orientation"))
            while(get_input(command, tag)) domain->erase_orientation(tag);
        else if(is_equal(object_type, "recorder"))
            while(get_input(command, tag)) domain->erase_recorder(tag);
        else if(is_equal(object_type, "section"))
            while(get_input(command, tag)) domain->erase_section(tag);
        else if(is_equal(object_type, "solver"))
            while(get_input(command, tag)) domain->erase_solver(tag);
        else if(is_equal(object_type, "step"))
            while(get_input(command, tag)) domain->erase_step(tag);

        return SUANPAN_SUCCESS;
    }

    int list_object(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        if(nullptr == domain) {
            suanpan_error("A valid domain is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string object_type;
        if(!get_input(command, object_type)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::vector<unsigned> list;
        if(is_equal(object_type, "amplitude"))
            for(auto& I : domain->get_amplitude_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "constraint"))
            for(auto& I : domain->get_constraint_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "converger"))
            for(auto& I : domain->get_converger_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "criterion"))
            for(auto& I : domain->get_criterion_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "element"))
            for(auto& I : domain->get_element_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "expression"))
            for(auto& I : domain->get_expression_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "group"))
            for(auto& I : domain->get_group_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "integrator"))
            for(auto& I : domain->get_integrator_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "load"))
            for(auto& I : domain->get_load_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "material"))
            for(auto& I : domain->get_material_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "modifier"))
            for(auto& I : domain->get_modifier_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "node"))
            for(auto& I : domain->get_node_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "orientation"))
            for(auto& I : domain->get_orientation_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "recorder"))
            for(auto& I : domain->get_recorder_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "section"))
            for(auto& I : domain->get_section_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "solver"))
            for(auto& I : domain->get_solver_pool()) list.emplace_back(I->get_tag());
        else if(is_equal(object_type, "step"))
            for(auto& I : domain->get_step_pool()) list.emplace_back(I.second->get_tag());

        suanpan_info("This domain has the following {}s:", object_type);
        for(const auto I : list)
            suanpan_info("\t{}", I);
        suanpan_info(".\n");

        return SUANPAN_SUCCESS;
    }

    int protect_object(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        if(nullptr == domain) {
            suanpan_error("A valid domain is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string object_type;
        if(!get_input(command, object_type)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(unsigned tag; is_equal(object_type, "element"))
            while(!command.eof() && get_input(command, tag)) {
                if(domain->find<Element>(tag)) domain->get<Element>(tag)->guard();
            }
        else if(is_equal(object_type, "node"))
            while(!command.eof() && get_input(command, tag)) {
                if(domain->find<Node>(tag)) domain->get<Node>(tag)->guard();
            }

        return SUANPAN_SUCCESS;
    }

    int save_object(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        if(nullptr == domain) {
            suanpan_error("A valid domain is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string object_id;
        if(!get_input(command, object_id)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(object_id, "Recorder")) {
            unsigned tag;
            while(get_input(command, tag))
                if(domain->find_recorder(tag)) domain->get_recorder(tag)->save();
        }
        else if(is_equal(object_id, "Stiffness")) {
            std::string name = "K";
            if(!command.eof() && !get_input(command, name)) name = "K";
            if(auto& stiffness = domain->get_factory()->get_stiffness()) stiffness->save(name.c_str());
        }
        else if(is_equal(object_id, "Mass")) {
            std::string name = "M";
            if(!command.eof() && !get_input(command, name)) name = "M";
            if(auto& mass = domain->get_factory()->get_mass()) mass->save(name.c_str());
        }
        else if(is_equal(object_id, "Damping")) {
            std::string name = "C";
            if(!command.eof() && !get_input(command, name)) name = "C";
            if(auto& damping = domain->get_factory()->get_damping()) damping->save(name.c_str());
        }
        else if(is_equal(object_id, "Model")) {
            std::string name = "Model.h5";
            if(!command.eof() && !get_input(command, name)) name = "Model.h5";
            domain->save(name);
        }

        return SUANPAN_SUCCESS;
    }

    int suspend_object(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        if(nullptr == domain) {
            suanpan_error("A valid domain is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string object_type;
        if(!get_input(command, object_type)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        const auto step_tag = domain->get_current_step_tag();

        if(unsigned tag; is_equal(object_type, "constraint"))
            while(!command.eof() && get_input(command, tag)) {
                if(domain->find_constraint(tag)) domain->get_constraint(tag)->set_end_step(step_tag);
            }
        else if(is_equal(object_type, "load"))
            while(!command.eof() && get_input(command, tag)) {
                if(domain->find_load(tag)) domain->get_load(tag)->set_end_step(step_tag);
            }

        return SUANPAN_SUCCESS;
    }

    int use_object(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        if(nullptr == domain) {
            suanpan_error("A valid domain is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string object_id;
        if(!get_input(command, object_id)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(object_id, "Integrator")) {
            unsigned tag;
            if(!get_input(command, tag)) {
                suanpan_error("A valid integrator tag is required.\n");
                return SUANPAN_SUCCESS;
            }
            if(domain->find_integrator(tag)) domain->set_current_integrator_tag(tag);
        }
        else if(is_equal(object_id, "Solver")) {
            unsigned tag;
            if(!get_input(command, tag)) {
                suanpan_error("A valid solver tag is required.\n");
                return SUANPAN_SUCCESS;
            }
            if(domain->find_solver(tag)) domain->set_current_solver_tag(tag);
        }
        else if(is_equal(object_id, "Converger")) {
            unsigned tag;
            if(!get_input(command, tag)) {
                suanpan_error("A valid converger tag is required.\n");
                return SUANPAN_SUCCESS;
            }
            if(domain->find_converger(tag)) domain->set_current_converger_tag(tag);
        }

        return SUANPAN_SUCCESS;
    }

    int set_property(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        std::string property_id;
        if(!get_input(command, property_id)) {
            suanpan_error("A valid property type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(property_id, "output_folder")) {
            std::string value;

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
        if(is_equal(property_id, "num_nodes")) {
            if(int value; get_input(command, value)) SUANPAN_NUM_NODES = value;
            else
                suanpan_error("A valid value is required.\n");

            return SUANPAN_SUCCESS;
        }
        if(is_equal(property_id, "screen_output")) {
            if(std::string value; get_input(command, value)) SUANPAN_PRINT = is_true(value);
            else
                suanpan_error("A valid value is required.\n");

            return SUANPAN_SUCCESS;
        }
        if(is_equal(property_id, "verbose_output")) {
            if(std::string value; get_input(command, value)) SUANPAN_VERBOSE = is_true(value);
            else
                suanpan_error("A valid value is required.\n");

            return SUANPAN_SUCCESS;
        }

        if(is_equal(property_id, "color_model")) {
            if(std::string value; !get_input(command, value))
                suanpan_error("A valid value is required.\n");
            else if(is_equal("WP", value)) domain->set_color_model(DomainBase::ColorMethod::WP);
            else if(is_equal("MIS", value)) domain->set_color_model(DomainBase::ColorMethod::MIS);
            else domain->set_color_model(DomainBase::ColorMethod::OFF);

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
            domain->set_attribute(DomainBase::ModalAttribute::LinearSystem);

            return SUANPAN_SUCCESS;
        }

        if(domain->get_current_step_tag() == 0) return SUANPAN_SUCCESS;

        auto& t_step = domain->get_current_step();

        if(is_equal(property_id, "fixed_step_size")) {
            std::string value;
            get_input(command, value) ? t_step->set_fixed_step_size(is_true(value)) : suanpan_error("A valid value is required.\n");
        }
        else if(is_equal(property_id, "symm_mat")) {
            std::string value;
            get_input(command, value) ? t_step->set_symm(is_true(value)) : suanpan_error("A valid value is required.\n");
        }
        else if(is_equal(property_id, "band_mat")) {
            std::string value;
            get_input(command, value) ? t_step->set_band(is_true(value)) : suanpan_error("A valid value is required.\n");
        }
        else if(is_equal(property_id, "sparse_mat")) {
            std::string value;
            get_input(command, value) ? t_step->set_sparse(is_true(value)) : suanpan_error("A valid value is required.\n");
        }
        else if(is_equal(property_id, "iterative_refinement")) {
            if(std::uint8_t value; get_input(command, value)) t_step->set_refinement(value);
            else
                suanpan_error("A valid value is required.\n");
        }
        else if(is_equal(property_id, "system_solver")) {
            if(std::string value; !get_input(command, value))
                suanpan_error("A valid value is required.\n");
            else if(is_equal(value, "LAPACK")) t_step->set_system_solver(SolverType::LAPACK);
            else if(is_equal(value, "SPIKE")) t_step->set_system_solver(SolverType::SPIKE);
            else if(is_equal(value, "SUPERLU")) t_step->set_system_solver(SolverType::SUPERLU);
            else if(is_equal(value, "LIS")) {
                t_step->set_system_solver(SolverType::LIS);
                t_step->set_lis_option(command);
            }
            else if(is_equal(value, "MUMPS")) {
                t_step->set_system_solver(SolverType::MUMPS);
                t_step->set_option(command);
            }
#ifdef SUANPAN_CUDA
            else if(is_equal(value, "CUDA")) t_step->set_system_solver(SolverType::CUDA);
#ifdef SUANPAN_MAGMA
            else if(is_equal(value, "MAGMA")) {
                t_step->set_system_solver(SolverType::MAGMA);
                t_step->set_option(command);
            }
#endif
#endif
#ifdef SUANPAN_MKL
            else if(is_equal(value, "PARDISO")) {
                t_step->set_system_solver(SolverType::PARDISO);
                t_step->set_option(command);
            }
            else if(is_equal(value, "FGMRES")) t_step->set_system_solver(SolverType::FGMRES);
#endif
            else
                suanpan_error("A valid solver type is required.\n");
        }
        else if(is_equal(property_id, "sub_system_solver")) {
            if(std::string value; !get_input(command, value))
                suanpan_error("A valid value is required.\n");
            else if(is_equal(value, "LAPACK")) t_step->set_sub_system_solver(SolverType::LAPACK);
            else if(is_equal(value, "SPIKE")) t_step->set_sub_system_solver(SolverType::SPIKE);
            else if(is_equal(value, "SUPERLU")) t_step->set_sub_system_solver(SolverType::SUPERLU);
            else if(is_equal(value, "LIS")) {
                t_step->set_sub_system_solver(SolverType::LIS);
                t_step->set_lis_option(command);
            }
            else if(is_equal(value, "MUMPS")) {
                t_step->set_sub_system_solver(SolverType::MUMPS);
                t_step->set_option(command);
            }
#ifdef SUANPAN_CUDA
            else if(is_equal(value, "CUDA")) t_step->set_sub_system_solver(SolverType::CUDA);
#endif
#ifdef SUANPAN_MKL
            else if(is_equal(value, "PARDISO")) {
                t_step->set_sub_system_solver(SolverType::PARDISO);
                t_step->set_option(command);
            }
            else if(is_equal(value, "FGMRES")) t_step->set_sub_system_solver(SolverType::FGMRES);
#endif
            else
                suanpan_error("A valid solver type is required.\n");
        }
        else if(is_equal(property_id, "precision")) {
            if(std::string value; !get_input(command, value))
                suanpan_error("A valid value is required.\n");
            else if(is_equal(value, "DOUBLE") || is_equal(value, "FULL")) t_step->set_precision(Precision::FULL);
            else if(is_equal(value, "SINGLE") || is_equal(value, "MIXED")) t_step->set_precision(Precision::MIXED);
            else
                suanpan_error("A valid precision is required.\n");
        }
        else if(is_equal(property_id, "tolerance") || is_equal(property_id, "fgmres_tolerance")) {
            double value;
            get_input(command, value) ? t_step->set_tolerance(value) : suanpan_error("A valid value is required.\n");
        }
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

    int print_info(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
        std::string object_type;
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
        else if(is_equal(object_type, "num_threads"))
            suanpan_info("SUANPAN_NUM_THREADS: {}\n", SUANPAN_NUM_THREADS);
        else if(is_equal(object_type, "num_nodes"))
            suanpan_info("SUANPAN_NUM_NODES: {}\n", SUANPAN_NUM_NODES);
        else if(is_equal(object_type, "statistics") || is_equal(object_type, "stats")) {
            suanpan_info("Updating element trial status used:\n\t{:.5E} s.\n", domain->stats<Statistics::UpdateStatus>());
            suanpan_info("Assembling global vector used:\n\t{:.5E} s.\n", domain->stats<Statistics::AssembleVector>());
            suanpan_info("Assembling global system used:\n\t{:.5E} s.\n", domain->stats<Statistics::AssembleMatrix>());
            suanpan_info("Processing constraints used:\n\t{:.5E} s.\n", domain->stats<Statistics::ProcessConstraint>());
            suanpan_info("Solving global system used:\n\t{:.5E} s.\n", domain->stats<Statistics::SolveSystem>());
        }

        return SUANPAN_SUCCESS;
    }

    int print_command() {
        // ReSharper disable once CppIfCanBeReplacedByConstexprIf
        // ReSharper disable once CppUnreachableCode
        if(0 != comm_rank) return SUANPAN_SUCCESS;

        suanpan_info("The available first-level commands are listed. Please check online manual for reference. https://tlcfem.github.io/suanPan-manual/latest/\n");

        constexpr auto format = "    {:>20}  {}\n";
        suanpan_info(format, "amplitude", "define amplitudes");
        suanpan_info(format, "analyze/analyse", "analyse the model");
        suanpan_info(format, "benchmark", "benchmark the platform for comparisons");
        suanpan_info(format, "clear", "clear the model to the initial state");
        suanpan_info(format, "command", "list all commands");
        suanpan_info(format, "constraint", "define constraints such as boundary conditions");
        suanpan_info(format, "converger", "define convergers");
        suanpan_info(format, "criterion", "define stopping criteria");
        suanpan_info(format, "delete/erase/remove", "delete objects");
        suanpan_info(format, "disable/mute", "disable objects");
        suanpan_info(format, "domain", "create/switch to other problem domains");
        suanpan_info(format, "element", "define elements");
        suanpan_info(format, "enable", "enable objects");
        suanpan_info(format, "example", "establish and execute a minimum example");
        suanpan_info(format, "exit/quit", "exit the program");
        suanpan_info(format, "expression", "define mathematical expressions to be used in other objects");
        suanpan_info(format, "file", "load external files");
        suanpan_info(format, "fullname", "print the full path of the program");
        suanpan_info(format, "group", "define groups via various rules");
        suanpan_info(format, "hdf5recorder", "define recorders using hdf5 format");
        suanpan_info(format, "import", "import external modules");
        suanpan_info(format, "initial", "define initial conditions for nodes and materials");
        suanpan_info(format, "integrator", "define time integration algorithms");
        suanpan_info(format, "list", "list objects in the current domain");
        suanpan_info(format, "load", "define loads of various types");
        suanpan_info(format, "material", "define materials");
        suanpan_info(format, "materialtest*", "test materials without creating a finite element model");
        suanpan_info(format, "modifier", "define modifiers that modify the existing model properties");
        suanpan_info(format, "node", "define nodes");
        suanpan_info(format, "orientation", "define beam section orientations");
        suanpan_info(format, "overview", "walk thorugh a quick overview of the application");
        suanpan_info(format, "peek", "peek the current information of the target object");
        suanpan_info(format, "plainrecorder", "define recorders using plain text format");
        suanpan_info(format, "plot", "plot and optionally save the model with VTK");
        suanpan_info(format, "precheck", "check the model without the actual analysis");
        suanpan_info(format, "protect", "protect objects from being disabled");
        suanpan_info(format, "pwd", "print/change the current working folder");
        suanpan_info(format, "qrcode", "print a qr code");
        suanpan_info(format, "recorder", "define recorders");
        suanpan_info(format, "reset", "reset the model to the previously converged state");
        suanpan_info(format, "response_spectrum", "compute the response spectrum of a given ground motion");
        suanpan_info(format, "save", "save objects");
        suanpan_info(format, "sdof_response", "compute the sdof response of a given ground motion");
        suanpan_info(format, "section", "define truss/beam/plate/shell sections");
        suanpan_info(format, "sectiontest*", "test sections without creating a finite element model");
        suanpan_info(format, "set", "set properties of the analysis/model");
        suanpan_info(format, "solver", "define solvers");
        suanpan_info(format, "step", "define analysis steps");
        suanpan_info(format, "summary", "print summary for the current problem domain");
        suanpan_info(format, "suspend", "suspend objects in the current step");
        suanpan_info(format, "terminal", "execute commands in terminal");
        suanpan_info(format, "upsampling", "upsample the given ground motion with various filters");
        suanpan_info(format, "version", "print version information");

        return SUANPAN_SUCCESS;
    }

    int run_example() {
        const auto new_model = std::make_shared<Bead>();

        suanpan_info("====================================================\n");
        suanpan_info("-> A Minimum Example: Elastic Truss Under Tension <-\n");
        suanpan_info("====================================================\n");

        constexpr auto wait_time = 1000;

        auto run_command = [&](const std::string& command_string) {
            suanpan_highlight("\t{}\n", command_string);
            auto command = std::istringstream(command_string);
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
} // namespace

int process_command(const shared_ptr<Bead>& model, std::istringstream& command) {
    if(nullptr == model) return SUANPAN_SUCCESS;

    std::string command_id;
    if(!get_input(command, command_id)) return SUANPAN_SUCCESS;

    if(is_equal(command_id, "exit") || is_equal(command_id, "quit")) return SUANPAN_EXIT;

    if(is_equal(command_id, "overview")) {
        overview();
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "file")) {
        std::string file_name;
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

    if(is_equal(command_id, "list")) return list_object(domain, command);
    if(is_equal(command_id, "protect")) return protect_object(domain, command);
    if(is_equal(command_id, "save")) return save_object(domain, command);
    if(is_equal(command_id, "set")) return set_property(domain, command);
    if(is_equal(command_id, "suspend")) return suspend_object(domain, command);
    if(is_equal(command_id, "use")) return use_object(domain, command);

    if(is_equal(command_id, "amplitude")) return create_new_amplitude(domain, command);
    if(is_equal(command_id, "expression")) return create_new_expression(domain, command);
    if(is_equal(command_id, "converger")) return create_new_converger(domain, command);
    if(is_equal(command_id, "criterion")) return create_new_criterion(domain, command);
    if(is_equal(command_id, "element")) return create_new_element(domain, command);
    if(is_equal(command_id, "hdf5recorder")) return create_new_recorder(domain, command, true);
    if(is_equal(command_id, "import")) return create_new_external_module(domain, command);
    if(is_equal(command_id, "initial")) return create_new_initial(domain, command);
    if(is_equal(command_id, "integrator")) return create_new_integrator(domain, command);
    if(is_equal(command_id, "mass")) return create_new_mass(domain, command);
    if(is_equal(command_id, "material")) return create_new_material(domain, command);
    if(is_equal(command_id, "modifier")) return create_new_modifier(domain, command);
    if(is_equal(command_id, "node")) return create_new_node(domain, command);
    if(is_equal(command_id, "orientation")) return create_new_orientation(domain, command);
    if(is_equal(command_id, "plainrecorder")) return create_new_recorder(domain, command, false);
    if(is_equal(command_id, "recorder")) return create_new_recorder(domain, command);
    if(is_equal(command_id, "section")) return create_new_section(domain, command);
    if(is_equal(command_id, "solver")) return create_new_solver(domain, command);
    if(is_equal(command_id, "step")) return create_new_step(domain, command);

    // groups
    auto group_handler = [&] {
        command.seekg(0);
        return create_new_group(domain, command);
    };

    if(is_equal(command_id, "group")) return create_new_group(domain, command);
    if(is_equal(command_id, "customnodegroup")) return group_handler();
    if(is_equal(command_id, "elementgroup")) return group_handler();
    if(is_equal(command_id, "generate")) return group_handler();
    if(is_equal(command_id, "generatebyplane")) return group_handler();
    if(is_equal(command_id, "generatebypoint")) return group_handler();
    if(is_equal(command_id, "generatebyrule")) return group_handler();
    if(is_equal(command_id, "groupgroup")) return group_handler();
    if(is_equal(command_id, "nodegroup")) return group_handler();

    // loads
    auto load_handler = [&] {
        command.seekg(0);
        return create_new_load(domain, command);
    };

    if(is_equal(command_id, "load")) return create_new_load(domain, command);
    if(is_equal(command_id, "acceleration")) return load_handler();
    if(is_equal(command_id, "bodyforce")) return load_handler();
    if(is_equal(command_id, "cload")) return load_handler();
    if(is_equal(command_id, "disp")) return load_handler();
    if(is_equal(command_id, "displacement")) return load_handler();
    if(is_equal(command_id, "dispload")) return load_handler();
    if(is_equal(command_id, "groupbodyforce")) return load_handler();
    if(is_equal(command_id, "groupcload")) return load_handler();
    if(is_equal(command_id, "groupdisp")) return load_handler();
    if(is_equal(command_id, "groupdisplacement")) return load_handler();
    if(is_equal(command_id, "groupdispload")) return load_handler();
    if(is_equal(command_id, "lineudl2d")) return load_handler();
    if(is_equal(command_id, "lineudl3d")) return load_handler();
    if(is_equal(command_id, "refereceload")) return load_handler();
    if(is_equal(command_id, "refforce")) return load_handler();
    if(is_equal(command_id, "refload")) return load_handler();
    if(is_equal(command_id, "supportacceleration")) return load_handler();
    if(is_equal(command_id, "supportdisplacement")) return load_handler();
    if(is_equal(command_id, "supportvelocity")) return load_handler();

    // constraints
    auto constraint_handler = [&] {
        command.seekg(0);
        return create_new_constraint(domain, command);
    };

    if(is_equal(command_id, "constraint")) return create_new_constraint(domain, command);
    if(is_equal(command_id, "embed2d") || is_equal(command_id, "embed3d")) return constraint_handler();
    if(is_equal(command_id, "finiterestitutionwall") || is_equal(command_id, "finiterestitutionwallpenalty")) return constraint_handler();
    if(is_equal(command_id, "finiterigidwall") || is_equal(command_id, "finiterigidwallpenalty")) return constraint_handler();
    if(is_equal(command_id, "finiterigidwallmultiplier")) return constraint_handler();
    if(is_equal(command_id, "fix")) return constraint_handler();
    if(is_equal(command_id, "fix2")) return constraint_handler();
    if(is_equal(command_id, "fixedlength2d")) return constraint_handler();
    if(is_equal(command_id, "fixedlength3d")) return constraint_handler();
    if(is_equal(command_id, "groupmultiplierbc")) return constraint_handler();
    if(is_equal(command_id, "grouppenaltybc")) return constraint_handler();
    if(is_equal(command_id, "maxforce2d") || is_equal(command_id, "maxforce3d")) return constraint_handler();
    if(is_equal(command_id, "maxgap2d") || is_equal(command_id, "maxgap3d")) return constraint_handler();
    if(is_equal(command_id, "mingap2d") || is_equal(command_id, "mingap3d")) return constraint_handler();
    if(is_equal(command_id, "mpc")) return constraint_handler();
    if(is_equal(command_id, "multiplierbc")) return constraint_handler();
    if(is_equal(command_id, "nodefacet")) return constraint_handler();
    if(is_equal(command_id, "nodeline")) return constraint_handler();
    if(is_equal(command_id, "particlecollision2d") || is_equal(command_id, "particlecollision3d")) return constraint_handler();
    if(is_equal(command_id, "penaltybc")) return constraint_handler();
    if(is_equal(command_id, "restitutionwall") || is_equal(command_id, "restitutionwallpenalty")) return constraint_handler();
    if(is_equal(command_id, "rigidwall") || is_equal(command_id, "rigidwallpenalty")) return constraint_handler();
    if(is_equal(command_id, "rigidwallmultiplier")) return constraint_handler();
    if(is_equal(command_id, "sleeve2d") || is_equal(command_id, "sleeve3d")) return constraint_handler();

    // testers
    if(is_equal(command_id, "materialtest1d")) return test_material(domain, command, 1);
    if(is_equal(command_id, "materialtest2d")) return test_material(domain, command, 3);
    if(is_equal(command_id, "materialtest3d")) return test_material(domain, command, 6);
    if(is_equal(command_id, "materialtestwithbase3d")) return test_material_with_base3d(domain, command);
    if(is_equal(command_id, "materialtestbyload1d")) return test_material_by_load(domain, command, 1);
    if(is_equal(command_id, "materialtestbyload2d")) return test_material_by_load(domain, command, 3);
    if(is_equal(command_id, "materialtestbyload3d")) return test_material_by_load(domain, command, 6);
    if(is_equal(command_id, "materialtestbyloadwithbase3d")) return test_material_by_load_with_base3d(domain, command);
    if(is_equal(command_id, "materialtestbystrainhistory")) return test_material_by_strain_history(domain, command);
    if(is_equal(command_id, "materialtestbystresshistory")) return test_material_by_stress_history(domain, command);

    if(is_equal(command_id, "sectiontest1d")) return test_section(domain, command, 1);
    if(is_equal(command_id, "sectiontest2d")) return test_section(domain, command, 2);
    if(is_equal(command_id, "sectiontest3d")) return test_section(domain, command, 3);
    if(is_equal(command_id, "sectiontestbydeformationhistory")) return test_section_by_deformation_history(domain, command);

    if(is_equal(command_id, "qrcode")) return qrcode();

    if(is_equal(command_id, "plot")) return vtk_parser(domain, command);

    if(is_equal(command_id, "peek")) return print_info(domain, command);

    if(is_equal(command_id, "command")) return print_command();

    if(is_equal(command_id, "example")) return run_example();

    if(is_equal(command_id, "precheck")) return model->precheck();

    if(is_equal(command_id, "analyze") || is_equal(command_id, "analyse")) {
        const auto options = get_remaining(command);
        if(SUANPAN_WARNING_COUNT > 0 && !if_contain(options, "ignore_warning") && !if_contain(options, "ignore-warning")) {
            suanpan_warning("There are {} warnings, please fix them first or use `ignore-warning` to ignore them.\n", SUANPAN_WARNING_COUNT);
            return SUANPAN_SUCCESS;
        }
        if(SUANPAN_ERROR_COUNT > 0 && !if_contain(options, "ignore_error") && !if_contain(options, "ignore-error")) {
            suanpan_warning("There are {} errors, please fix them first or use `ignore-error` to ignore them.\n", SUANPAN_ERROR_COUNT);
            return SUANPAN_SUCCESS;
        }
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
        else if(std::string path; get_input(command, path)) {
            std::error_code code;
            fs::current_path(path, code);
            if(0 != code.value())
                suanpan_error("Fail to set path \"{}\"\n", code.category().message(code.value()));
        }
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "sleep")) {
        if(auto t = 1000ll; get_optional_input(command, t)) std::this_thread::sleep_for(std::chrono::milliseconds(t));
        else
            suanpan_error("A positive integer in milliseconds is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "benchmark")) return benchmark();

    if(is_equal(command_id, "clear")) {
        domain->wait();

        auto flag = true;
        for(auto& t_integrator : domain->get_integrator_pool())
            if(nullptr != t_integrator->get_domain()) {
                t_integrator->clear_status();
                flag = false;
            }

        if(flag) domain->clear_status();

        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "reset")) {
        domain->wait();

        auto flag = true;
        for(auto& t_integrator : domain->get_integrator_pool())
            if(nullptr != t_integrator->get_domain()) {
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

    if(is_equal(command_id, "terminal") || is_equal(command_id, "t")) {
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

    if(is_equal(command_id, "license")) {
        suanpan_info(R"(                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

                            Preamble

  The GNU General Public License is a free, copyleft license for
software and other kinds of works.

  The licenses for most software and other practical works are designed
to take away your freedom to share and change the works.  By contrast,
the GNU General Public License is intended to guarantee your freedom to
share and change all versions of a program--to make sure it remains free
software for all its users.  We, the Free Software Foundation, use the
GNU General Public License for most of our software; it applies also to
any other work released this way by its authors.  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
them if you wish), that you receive source code or can get it if you
want it, that you can change the software or use pieces of it in new
free programs, and that you know you can do these things.

  To protect your rights, we need to prevent others from denying you
these rights or asking you to surrender the rights.  Therefore, you have
certain responsibilities if you distribute copies of the software, or if
you modify it: responsibilities to respect the freedom of others.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must pass on to the recipients the same
freedoms that you received.  You must make sure that they, too, receive
or can get the source code.  And you must show them these terms so they
know their rights.

  Developers that use the GNU GPL protect your rights with two steps:
(1) assert copyright on the software, and (2) offer you this License
giving you legal permission to copy, distribute and/or modify it.

  For the developers' and authors' protection, the GPL clearly explains
that there is no warranty for this free software.  For both users' and
authors' sake, the GPL requires that modified versions be marked as
changed, so that their problems will not be attributed erroneously to
authors of previous versions.

  Some devices are designed to deny users access to install or run
modified versions of the software inside them, although the manufacturer
can do so.  This is fundamentally incompatible with the aim of
protecting users' freedom to change the software.  The systematic
pattern of such abuse occurs in the area of products for individuals to
use, which is precisely where it is most unacceptable.  Therefore, we
have designed this version of the GPL to prohibit the practice for those
products.  If such problems arise substantially in other domains, we
stand ready to extend this provision to those domains in future versions
of the GPL, as needed to protect the freedom of users.

  Finally, every program is threatened constantly by software patents.
States should not allow patents to restrict development and use of
software on general-purpose computers, but in those that do, we wish to
avoid the special danger that patents applied to a free program could
make it effectively proprietary.  To prevent this, the GPL assures that
patents cannot be used to render the program non-free.

  The precise terms and conditions for copying, distribution and
modification follow.

                       TERMS AND CONDITIONS

  0. Definitions.

  "This License" refers to version 3 of the GNU General Public License.

  "Copyright" also means copyright-like laws that apply to other kinds of
works, such as semiconductor masks.

  "The Program" refers to any copyrightable work licensed under this
License.  Each licensee is addressed as "you".  "Licensees" and
"recipients" may be individuals or organizations.

  To "modify" a work means to copy from or adapt all or part of the work
in a fashion requiring copyright permission, other than the making of an
exact copy.  The resulting work is called a "modified version" of the
earlier work or a work "based on" the earlier work.

  A "covered work" means either the unmodified Program or a work based
on the Program.

  To "propagate" a work means to do anything with it that, without
permission, would make you directly or secondarily liable for
infringement under applicable copyright law, except executing it on a
computer or modifying a private copy.  Propagation includes copying,
distribution (with or without modification), making available to the
public, and in some countries other activities as well.

  To "convey" a work means any kind of propagation that enables other
parties to make or receive copies.  Mere interaction with a user through
a computer network, with no transfer of a copy, is not conveying.

  An interactive user interface displays "Appropriate Legal Notices"
to the extent that it includes a convenient and prominently visible
feature that (1) displays an appropriate copyright notice, and (2)
tells the user that there is no warranty for the work (except to the
extent that warranties are provided), that licensees may convey the
work under this License, and how to view a copy of this License.  If
the interface presents a list of user commands or options, such as a
menu, a prominent item in the list meets this criterion.

  1. Source Code.

  The "source code" for a work means the preferred form of the work
for making modifications to it.  "Object code" means any non-source
form of a work.

  A "Standard Interface" means an interface that either is an official
standard defined by a recognized standards body, or, in the case of
interfaces specified for a particular programming language, one that
is widely used among developers working in that language.

  The "System Libraries" of an executable work include anything, other
than the work as a whole, that (a) is included in the normal form of
packaging a Major Component, but which is not part of that Major
Component, and (b) serves only to enable use of the work with that
Major Component, or to implement a Standard Interface for which an
implementation is available to the public in source code form.  A
"Major Component", in this context, means a major essential component
(kernel, window system, and so on) of the specific operating system
(if any) on which the executable work runs, or a compiler used to
produce the work, or an object code interpreter used to run it.

  The "Corresponding Source" for a work in object code form means all
the source code needed to generate, install, and (for an executable
work) run the object code and to modify the work, including scripts to
control those activities.  However, it does not include the work's
System Libraries, or general-purpose tools or generally available free
programs which are used unmodified in performing those activities but
which are not part of the work.  For example, Corresponding Source
includes interface definition files associated with source files for
the work, and the source code for shared libraries and dynamically
linked subprograms that the work is specifically designed to require,
such as by intimate data communication or control flow between those
subprograms and other parts of the work.

  The Corresponding Source need not include anything that users
can regenerate automatically from other parts of the Corresponding
Source.

  The Corresponding Source for a work in source code form is that
same work.

  2. Basic Permissions.

  All rights granted under this License are granted for the term of
copyright on the Program, and are irrevocable provided the stated
conditions are met.  This License explicitly affirms your unlimited
permission to run the unmodified Program.  The output from running a
covered work is covered by this License only if the output, given its
content, constitutes a covered work.  This License acknowledges your
rights of fair use or other equivalent, as provided by copyright law.

  You may make, run and propagate covered works that you do not
convey, without conditions so long as your license otherwise remains
in force.  You may convey covered works to others for the sole purpose
of having them make modifications exclusively for you, or provide you
with facilities for running those works, provided that you comply with
the terms of this License in conveying all material for which you do
not control copyright.  Those thus making or running the covered works
for you must do so exclusively on your behalf, under your direction
and control, on terms that prohibit them from making any copies of
your copyrighted material outside their relationship with you.

  Conveying under any other circumstances is permitted solely under
the conditions stated below.  Sublicensing is not allowed; section 10
makes it unnecessary.

  3. Protecting Users' Legal Rights From Anti-Circumvention Law.

  No covered work shall be deemed part of an effective technological
measure under any applicable law fulfilling obligations under article
11 of the WIPO copyright treaty adopted on 20 December 1996, or
similar laws prohibiting or restricting circumvention of such
measures.

  When you convey a covered work, you waive any legal power to forbid
circumvention of technological measures to the extent such circumvention
is effected by exercising rights under this License with respect to
the covered work, and you disclaim any intention to limit operation or
modification of the work as a means of enforcing, against the work's
users, your or third parties' legal rights to forbid circumvention of
technological measures.

  4. Conveying Verbatim Copies.

  You may convey verbatim copies of the Program's source code as you
receive it, in any medium, provided that you conspicuously and
appropriately publish on each copy an appropriate copyright notice;
keep intact all notices stating that this License and any
non-permissive terms added in accord with section 7 apply to the code;
keep intact all notices of the absence of any warranty; and give all
recipients a copy of this License along with the Program.

  You may charge any price or no price for each copy that you convey,
and you may offer support or warranty protection for a fee.

  5. Conveying Modified Source Versions.

  You may convey a work based on the Program, or the modifications to
produce it from the Program, in the form of source code under the
terms of section 4, provided that you also meet all of these conditions:

    a) The work must carry prominent notices stating that you modified
    it, and giving a relevant date.

    b) The work must carry prominent notices stating that it is
    released under this License and any conditions added under section
    7.  This requirement modifies the requirement in section 4 to
    "keep intact all notices".

    c) You must license the entire work, as a whole, under this
    License to anyone who comes into possession of a copy.  This
    License will therefore apply, along with any applicable section 7
    additional terms, to the whole of the work, and all its parts,
    regardless of how they are packaged.  This License gives no
    permission to license the work in any other way, but it does not
    invalidate such permission if you have separately received it.

    d) If the work has interactive user interfaces, each must display
    Appropriate Legal Notices; however, if the Program has interactive
    interfaces that do not display Appropriate Legal Notices, your
    work need not make them do so.

  A compilation of a covered work with other separate and independent
works, which are not by their nature extensions of the covered work,
and which are not combined with it such as to form a larger program,
in or on a volume of a storage or distribution medium, is called an
"aggregate" if the compilation and its resulting copyright are not
used to limit the access or legal rights of the compilation's users
beyond what the individual works permit.  Inclusion of a covered work
in an aggregate does not cause this License to apply to the other
parts of the aggregate.

  6. Conveying Non-Source Forms.

  You may convey a covered work in object code form under the terms
of sections 4 and 5, provided that you also convey the
machine-readable Corresponding Source under the terms of this License,
in one of these ways:

    a) Convey the object code in, or embodied in, a physical product
    (including a physical distribution medium), accompanied by the
    Corresponding Source fixed on a durable physical medium
    customarily used for software interchange.

    b) Convey the object code in, or embodied in, a physical product
    (including a physical distribution medium), accompanied by a
    written offer, valid for at least three years and valid for as
    long as you offer spare parts or customer support for that product
    model, to give anyone who possesses the object code either (1) a
    copy of the Corresponding Source for all the software in the
    product that is covered by this License, on a durable physical
    medium customarily used for software interchange, for a price no
    more than your reasonable cost of physically performing this
    conveying of source, or (2) access to copy the
    Corresponding Source from a network server at no charge.

    c) Convey individual copies of the object code with a copy of the
    written offer to provide the Corresponding Source.  This
    alternative is allowed only occasionally and noncommercially, and
    only if you received the object code with such an offer, in accord
    with subsection 6b.

    d) Convey the object code by offering access from a designated
    place (gratis or for a charge), and offer equivalent access to the
    Corresponding Source in the same way through the same place at no
    further charge.  You need not require recipients to copy the
    Corresponding Source along with the object code.  If the place to
    copy the object code is a network server, the Corresponding Source
    may be on a different server (operated by you or a third party)
    that supports equivalent copying facilities, provided you maintain
    clear directions next to the object code saying where to find the
    Corresponding Source.  Regardless of what server hosts the
    Corresponding Source, you remain obligated to ensure that it is
    available for as long as needed to satisfy these requirements.

    e) Convey the object code using peer-to-peer transmission, provided
    you inform other peers where the object code and Corresponding
    Source of the work are being offered to the general public at no
    charge under subsection 6d.

  A separable portion of the object code, whose source code is excluded
from the Corresponding Source as a System Library, need not be
included in conveying the object code work.

  A "User Product" is either (1) a "consumer product", which means any
tangible personal property which is normally used for personal, family,
or household purposes, or (2) anything designed or sold for incorporation
into a dwelling.  In determining whether a product is a consumer product,
doubtful cases shall be resolved in favor of coverage.  For a particular
product received by a particular user, "normally used" refers to a
typical or common use of that class of product, regardless of the status
of the particular user or of the way in which the particular user
actually uses, or expects or is expected to use, the product.  A product
is a consumer product regardless of whether the product has substantial
commercial, industrial or non-consumer uses, unless such uses represent
the only significant mode of use of the product.

  "Installation Information" for a User Product means any methods,
procedures, authorization keys, or other information required to install
and execute modified versions of a covered work in that User Product from
a modified version of its Corresponding Source.  The information must
suffice to ensure that the continued functioning of the modified object
code is in no case prevented or interfered with solely because
modification has been made.)"
                     R"(

  If you convey an object code work under this section in, or with, or
specifically for use in, a User Product, and the conveying occurs as
part of a transaction in which the right of possession and use of the
User Product is transferred to the recipient in perpetuity or for a
fixed term (regardless of how the transaction is characterized), the
Corresponding Source conveyed under this section must be accompanied
by the Installation Information.  But this requirement does not apply
if neither you nor any third party retains the ability to install
modified object code on the User Product (for example, the work has
been installed in ROM).

  The requirement to provide Installation Information does not include a
requirement to continue to provide support service, warranty, or updates
for a work that has been modified or installed by the recipient, or for
the User Product in which it has been modified or installed.  Access to a
network may be denied when the modification itself materially and
adversely affects the operation of the network or violates the rules and
protocols for communication across the network.

  Corresponding Source conveyed, and Installation Information provided,
in accord with this section must be in a format that is publicly
documented (and with an implementation available to the public in
source code form), and must require no special password or key for
unpacking, reading or copying.

  7. Additional Terms.

  "Additional permissions" are terms that supplement the terms of this
License by making exceptions from one or more of its conditions.
Additional permissions that are applicable to the entire Program shall
be treated as though they were included in this License, to the extent
that they are valid under applicable law.  If additional permissions
apply only to part of the Program, that part may be used separately
under those permissions, but the entire Program remains governed by
this License without regard to the additional permissions.

  When you convey a copy of a covered work, you may at your option
remove any additional permissions from that copy, or from any part of
it.  (Additional permissions may be written to require their own
removal in certain cases when you modify the work.)  You may place
additional permissions on material, added by you to a covered work,
for which you have or can give appropriate copyright permission.

  Notwithstanding any other provision of this License, for material you
add to a covered work, you may (if authorized by the copyright holders of
that material) supplement the terms of this License with terms:

    a) Disclaiming warranty or limiting liability differently from the
    terms of sections 15 and 16 of this License; or

    b) Requiring preservation of specified reasonable legal notices or
    author attributions in that material or in the Appropriate Legal
    Notices displayed by works containing it; or

    c) Prohibiting misrepresentation of the origin of that material, or
    requiring that modified versions of such material be marked in
    reasonable ways as different from the original version; or

    d) Limiting the use for publicity purposes of names of licensors or
    authors of the material; or

    e) Declining to grant rights under trademark law for use of some
    trade names, trademarks, or service marks; or

    f) Requiring indemnification of licensors and authors of that
    material by anyone who conveys the material (or modified versions of
    it) with contractual assumptions of liability to the recipient, for
    any liability that these contractual assumptions directly impose on
    those licensors and authors.

  All other non-permissive additional terms are considered "further
restrictions" within the meaning of section 10.  If the Program as you
received it, or any part of it, contains a notice stating that it is
governed by this License along with a term that is a further
restriction, you may remove that term.  If a license document contains
a further restriction but permits relicensing or conveying under this
License, you may add to a covered work material governed by the terms
of that license document, provided that the further restriction does
not survive such relicensing or conveying.

  If you add terms to a covered work in accord with this section, you
must place, in the relevant source files, a statement of the
additional terms that apply to those files, or a notice indicating
where to find the applicable terms.

  Additional terms, permissive or non-permissive, may be stated in the
form of a separately written license, or stated as exceptions;
the above requirements apply either way.

  8. Termination.

  You may not propagate or modify a covered work except as expressly
provided under this License.  Any attempt otherwise to propagate or
modify it is void, and will automatically terminate your rights under
this License (including any patent licenses granted under the third
paragraph of section 11).

  However, if you cease all violation of this License, then your
license from a particular copyright holder is reinstated (a)
provisionally, unless and until the copyright holder explicitly and
finally terminates your license, and (b) permanently, if the copyright
holder fails to notify you of the violation by some reasonable means
prior to 60 days after the cessation.

  Moreover, your license from a particular copyright holder is
reinstated permanently if the copyright holder notifies you of the
violation by some reasonable means, this is the first time you have
received notice of violation of this License (for any work) from that
copyright holder, and you cure the violation prior to 30 days after
your receipt of the notice.

  Termination of your rights under this section does not terminate the
licenses of parties who have received copies or rights from you under
this License.  If your rights have been terminated and not permanently
reinstated, you do not qualify to receive new licenses for the same
material under section 10.

  9. Acceptance Not Required for Having Copies.

  You are not required to accept this License in order to receive or
run a copy of the Program.  Ancillary propagation of a covered work
occurring solely as a consequence of using peer-to-peer transmission
to receive a copy likewise does not require acceptance.  However,
nothing other than this License grants you permission to propagate or
modify any covered work.  These actions infringe copyright if you do
not accept this License.  Therefore, by modifying or propagating a
covered work, you indicate your acceptance of this License to do so.

  10. Automatic Licensing of Downstream Recipients.

  Each time you convey a covered work, the recipient automatically
receives a license from the original licensors, to run, modify and
propagate that work, subject to this License.  You are not responsible
for enforcing compliance by third parties with this License.

  An "entity transaction" is a transaction transferring control of an
organization, or substantially all assets of one, or subdividing an
organization, or merging organizations.  If propagation of a covered
work results from an entity transaction, each party to that
transaction who receives a copy of the work also receives whatever
licenses to the work the party's predecessor in interest had or could
give under the previous paragraph, plus a right to possession of the
Corresponding Source of the work from the predecessor in interest, if
the predecessor has it or can get it with reasonable efforts.

  You may not impose any further restrictions on the exercise of the
rights granted or affirmed under this License.  For example, you may
not impose a license fee, royalty, or other charge for exercise of
rights granted under this License, and you may not initiate litigation
(including a cross-claim or counterclaim in a lawsuit) alleging that
any patent claim is infringed by making, using, selling, offering for
sale, or importing the Program or any portion of it.

  11. Patents.

  A "contributor" is a copyright holder who authorizes use under this
License of the Program or a work on which the Program is based.  The
work thus licensed is called the contributor's "contributor version".

  A contributor's "essential patent claims" are all patent claims
owned or controlled by the contributor, whether already acquired or
hereafter acquired, that would be infringed by some manner, permitted
by this License, of making, using, or selling its contributor version,
but do not include claims that would be infringed only as a
consequence of further modification of the contributor version.  For
purposes of this definition, "control" includes the right to grant
patent sublicenses in a manner consistent with the requirements of
this License.

  Each contributor grants you a non-exclusive, worldwide, royalty-free
patent license under the contributor's essential patent claims, to
make, use, sell, offer for sale, import and otherwise run, modify and
propagate the contents of its contributor version.

  In the following three paragraphs, a "patent license" is any express
agreement or commitment, however denominated, not to enforce a patent
(such as an express permission to practice a patent or covenant not to
sue for patent infringement).  To "grant" such a patent license to a
party means to make such an agreement or commitment not to enforce a
patent against the party.

  If you convey a covered work, knowingly relying on a patent license,
and the Corresponding Source of the work is not available for anyone
to copy, free of charge and under the terms of this License, through a
publicly available network server or other readily accessible means,
then you must either (1) cause the Corresponding Source to be so
available, or (2) arrange to deprive yourself of the benefit of the
patent license for this particular work, or (3) arrange, in a manner
consistent with the requirements of this License, to extend the patent
license to downstream recipients.  "Knowingly relying" means you have
actual knowledge that, but for the patent license, your conveying the
covered work in a country, or your recipient's use of the covered work
in a country, would infringe one or more identifiable patents in that
country that you have reason to believe are valid.

  If, pursuant to or in connection with a single transaction or
arrangement, you convey, or propagate by procuring conveyance of, a
covered work, and grant a patent license to some of the parties
receiving the covered work authorizing them to use, propagate, modify
or convey a specific copy of the covered work, then the patent license
you grant is automatically extended to all recipients of the covered
work and works based on it.

  A patent license is "discriminatory" if it does not include within
the scope of its coverage, prohibits the exercise of, or is
conditioned on the non-exercise of one or more of the rights that are
specifically granted under this License.  You may not convey a covered
work if you are a party to an arrangement with a third party that is
in the business of distributing software, under which you make payment
to the third party based on the extent of your activity of conveying
the work, and under which the third party grants, to any of the
parties who would receive the covered work from you, a discriminatory
patent license (a) in connection with copies of the covered work
conveyed by you (or copies made from those copies), or (b) primarily
for and in connection with specific products or compilations that
contain the covered work, unless you entered into that arrangement,
or that patent license was granted, prior to 28 March 2007.

  Nothing in this License shall be construed as excluding or limiting
any implied license or other defenses to infringement that may
otherwise be available to you under applicable patent law.

  12. No Surrender of Others' Freedom.

  If conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot convey a
covered work so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you may
not convey it at all.  For example, if you agree to terms that obligate you
to collect a royalty for further conveying from those to whom you convey
the Program, the only way you could satisfy both those terms and this
License would be to refrain entirely from conveying the Program.

  13. Use with the GNU Affero General Public License.

  Notwithstanding any other provision of this License, you have
permission to link or combine any covered work with a work licensed
under version 3 of the GNU Affero General Public License into a single
combined work, and to convey the resulting work.  The terms of this
License will continue to apply to the part which is the covered work,
but the special requirements of the GNU Affero General Public License,
section 13, concerning interaction through a network will apply to the
combination as such.

  14. Revised Versions of this License.

  The Free Software Foundation may publish revised and/or new versions of
the GNU General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

  Each version is given a distinguishing version number.  If the
Program specifies that a certain numbered version of the GNU General
Public License "or any later version" applies to it, you have the
option of following the terms and conditions either of that numbered
version or of any later version published by the Free Software
Foundation.  If the Program does not specify a version number of the
GNU General Public License, you may choose any version ever published
by the Free Software Foundation.

  If the Program specifies that a proxy can decide which future
versions of the GNU General Public License can be used, that proxy's
public statement of acceptance of a version permanently authorizes you
to choose that version for the Program.

  Later license versions may give you additional or different
permissions.  However, no additional obligations are imposed on any
author or copyright holder as a result of your choosing to follow a
later version.

  15. Disclaimer of Warranty.

  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

  16. Limitation of Liability.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

  17. Interpretation of Sections 15 and 16.

  If the disclaimer of warranty and limitation of liability provided
above cannot be given local legal effect according to their terms,
reviewing courts shall apply local law that most closely approximates
an absolute waiver of all civil liability in connection with the
Program, unless a warranty or assumption of liability accompanies a
copy of the Program in return for a fee.
)");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "version")) print_version();
    else
        suanpan_error("Command \"{}\" not found.\n", command.str());

    return SUANPAN_SUCCESS;
}

auto sanitize_command(std::string& line, const char delimiter) {
    while(!line.empty() && delimiter == line.back()) line.pop_back();
}

bool normalise_command(std::string& all_line, std::string& command_line) {
    // if to parse and process immediately
    auto process = true;

    // clear comment line
    if(!command_line.empty() && '#' == command_line.front()) command_line.clear();
    // remove inline comment
    if(const auto if_comment = command_line.find('!'); std::string::npos != if_comment) command_line.erase(if_comment);
    // remove all delimiters
    for(auto& c : command_line)
        if(',' == c || '\t' == c || '\r' == c || '\n' == c) c = ' ';
    sanitize_command(command_line, ' ');
    // it is a command spanning multiple lines
    if(!command_line.empty() && '\\' == command_line.back()) {
        sanitize_command(command_line, '\\');
        process = false;
    }
    all_line.append(command_line);

    if(process) sanitize_command(all_line, ' ');

    return process && !all_line.empty();
}

int process_file(const shared_ptr<Bead>& model, const char* file_name) {
    std::vector<std::string> file_list;
    file_list.reserve(9);

    std::string str_name(file_name);
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

    uintmax_t input_file_size{};

    for(const auto& file : file_list) {
        const auto file_path = fs::path(file);
        input_file.open(file_path);
        if(input_file.is_open()) {
            input_file_size = file_size(file_path);
            break;
        }
    }

    if(!input_file.is_open()) {
        suanpan_error("Cannot open the input file \"{}\".\n", fs::path(file_name).generic_string());
        return SUANPAN_EXIT;
    }

    const auto history_path = get_history_path();
    ofstream output_file(history_path, !exists(history_path) || file_size(history_path) > 16777216 ? std::ios_base::out : (std::ios_base::app | std::ios_base::out));

    const auto record_command = output_file.is_open() && input_file_size <= 102400;

    if(record_command) output_file << "### start processing --> " << file_name << '\n';

    std::string all_line, command_line;
    while(!getline(input_file, command_line).fail()) {
        if(!normalise_command(all_line, command_line)) continue;
        // now process the command
        if(record_command) output_file << all_line << '\n';
        if(std::istringstream tmp_str(all_line); process_command(model, tmp_str) == SUANPAN_EXIT) {
            if(record_command) output_file << "### finish processing --> " << file_name << '\n';
            return SUANPAN_EXIT;
        }
        all_line.clear();
    }

    if(record_command) output_file << "### finish processing --> " << file_name << '\n';
    return SUANPAN_SUCCESS;
}

int execute_command(std::istringstream& command) {
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

#ifdef SUANPAN_MSVC
#pragma warning(disable : 4996)
#endif

fs::path get_history_path() {
#ifdef SUANPAN_WIN
    // ReSharper disable once CppDeprecatedEntity
    auto history_path = fs::path(getenv("USERPROFILE")); // NOLINT(concurrency-mt-unsafe, clang-diagnostic-deprecated-declarations)
#else
    auto history_path = fs::path(getenv("HOME"));
#endif

    history_path.append(".suanpan-history.sp");

    return history_path;
}
