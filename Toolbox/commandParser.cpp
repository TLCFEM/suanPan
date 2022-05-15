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

// ReSharper disable CppClangTidyCppcoreguidelinesInitVariables
#include <suanPan>
#include <thread>
#ifdef SUANPAN_WIN
#include <Windows.h>
#endif

using std::ifstream;
using std::string;
using std::vector;

int SUANPAN_NUM_THREADS = 10;
fs::path SUANPAN_OUTPUT = fs::current_path();

int process_command(const shared_ptr<Bead>& model, istringstream& command) {
    if(nullptr == model) return SUANPAN_SUCCESS;

    string command_id;
    if(!get_input(command, command_id)) return SUANPAN_SUCCESS;

    if(is_equal(command_id, "exit")) return SUANPAN_EXIT;
    if(is_equal(command_id, "quit")) return SUANPAN_EXIT;

    if(is_equal(command_id, "file")) return process_file(model, command);
    if(is_equal(command_id, "load")) return process_file(model, command);

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
    if(is_equal(command_id, "acceleration")) return create_new_acceleration(domain, command);
    if(is_equal(command_id, "amplitude")) return create_new_amplitude(domain, command);
    if(is_equal(command_id, "bodyforce")) return create_new_bodyforce(domain, command, false);
    if(is_equal(command_id, "cload")) return create_new_cload(domain, command);
    if(is_equal(command_id, "lineudl2d")) return create_new_lineudl(domain, command, 2);
    if(is_equal(command_id, "lineudl3d")) return create_new_lineudl(domain, command, 3);
    if(is_equal(command_id, "converger")) return create_new_converger(domain, command);
    if(is_equal(command_id, "constraint")) return create_new_constraint(domain, command);
    if(is_equal(command_id, "criterion")) return create_new_criterion(domain, command);
    if(is_equal(command_id, "disp")) return create_new_displacement(domain, command);
    if(is_equal(command_id, "displacement")) return create_new_displacement(domain, command);
    if(is_equal(command_id, "dispload")) return create_new_displacement(domain, command);
    if(is_equal(command_id, "element")) return create_new_element(domain, command);
    if(is_equal(command_id, "elementgroup")) return create_new_elementgroup(domain, command);
    if(is_equal(command_id, "finiterigidwall")) return create_new_rigidwall(domain, command, true, true);
    if(is_equal(command_id, "finiterigidwallmultiplier")) return create_new_rigidwall(domain, command, true, false);
    if(is_equal(command_id, "fix")) return create_new_bc(domain, command, true);
    if(is_equal(command_id, "fix2")) return create_new_bc(domain, command, false);
    if(is_equal(command_id, "fixedlength2d")) return create_new_fixedlength(domain, command, 2);
    if(is_equal(command_id, "fixedlength3d")) return create_new_fixedlength(domain, command, 3);
    if(is_equal(command_id, "generate")) return create_new_generate(domain, command);
    if(is_equal(command_id, "generatebyrule")) return create_new_generatebyrule(domain, command);
    if(is_equal(command_id, "generatebypoint")) return create_new_generatebypoint(domain, command);
    if(is_equal(command_id, "generatebyplane")) return create_new_generatebyplane(domain, command);
    if(is_equal(command_id, "groupbodyforce")) return create_new_bodyforce(domain, command, true);
    if(is_equal(command_id, "groupcload")) return create_new_cload(domain, command, true);
    if(is_equal(command_id, "groupdisp")) return create_new_displacement(domain, command, true);
    if(is_equal(command_id, "groupdisplacement")) return create_new_displacement(domain, command, true);
    if(is_equal(command_id, "groupdispload")) return create_new_displacement(domain, command, true);
    if(is_equal(command_id, "groupgroup")) return create_new_groupgroup(domain, command);
    if(is_equal(command_id, "groupmultiplierbc")) return create_new_groupbc(domain, command, false);
    if(is_equal(command_id, "grouppenaltybc")) return create_new_groupbc(domain, command, true);
    if(is_equal(command_id, "hdf5recorder")) return create_new_hdf5recorder(domain, command);
    if(is_equal(command_id, "import")) return create_new_external_module(domain, command);
    if(is_equal(command_id, "initial")) return create_new_initial(domain, command);
    if(is_equal(command_id, "integrator")) return create_new_integrator(domain, command);
    if(is_equal(command_id, "mass")) return create_new_mass(domain, command);
    if(is_equal(command_id, "material")) return create_new_material(domain, command);
    if(is_equal(command_id, "modifier")) return create_new_modifier(domain, command);
    if(is_equal(command_id, "mpc")) return create_new_mpc(domain, command);
    if(is_equal(command_id, "multiplierbc")) return create_new_bc(domain, command, false);
    if(is_equal(command_id, "node")) return create_new_node(domain, command);
    if(is_equal(command_id, "nodegroup")) return create_new_nodegroup(domain, command);
    if(is_equal(command_id, "orientation")) return create_new_orientation(domain, command);
    if(is_equal(command_id, "particlecollision2d")) return create_new_particlecollision2d(domain, command);
    if(is_equal(command_id, "particlecollision3d")) return create_new_particlecollision3d(domain, command);
    if(is_equal(command_id, "penaltybc")) return create_new_bc(domain, command, true);
    if(is_equal(command_id, "plainrecorder")) return create_new_plainrecorder(domain, command);
    if(is_equal(command_id, "recorder")) return create_new_recorder(domain, command);
    if(is_equal(command_id, "rigidwall")) return create_new_rigidwall(domain, command, false, true);
    if(is_equal(command_id, "rigidwallmultiplier")) return create_new_rigidwall(domain, command, false, false);
    if(is_equal(command_id, "section")) return create_new_section(domain, command);
    if(is_equal(command_id, "solver")) return create_new_solver(domain, command);
    if(is_equal(command_id, "step")) return create_new_step(domain, command);
    if(is_equal(command_id, "supportdisplacement")) return create_new_supportmotion(domain, command, 0);
    if(is_equal(command_id, "supportvelocity")) return create_new_supportmotion(domain, command, 1);
    if(is_equal(command_id, "supportacceleration")) return create_new_supportmotion(domain, command, 2);
    if(is_equal(command_id, "set")) return set_property(domain, command);

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

    if(is_equal(command_id, "analyze")) {
        const auto code = model->analyze();
        suanpan_info("\n");
        return code;
    }
    if(is_equal(command_id, "analyse")) {
        const auto code = model->analyze();
        suanpan_info("\n");
        return code;
    }

    if(is_equal(command_id, "fullname")) {
        suanpan_info("%s\n", SUANPAN_EXE);
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "pwd")) {
        if(command.eof()) suanpan_info("%s\n", fs::current_path().generic_string().c_str());
        else if(string path; get_input(command, path)) {
            std::error_code code;
            fs::current_path(path, code);
            if(0 != code.value()) suanpan_error("fail to set current path: %s\n", code.category().message(code.value()).c_str());
        }
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "help")) {
        print_helper();
        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "clear")) {
        domain->wait();

        auto flag = true;
        for(const auto& t_integrator : domain->get_integrator_pool())
            if(t_integrator->get_domain().lock() != nullptr) {
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
            if(t_integrator->get_domain().lock() != nullptr) {
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
#ifdef SUANPAN_WIN
        const auto handle = GetStdHandle(STD_OUTPUT_HANDLE);
        CONSOLE_SCREEN_BUFFER_INFO info;
        GetConsoleScreenBufferInfo(handle, &info);
        const auto current_attribute = info.wAttributes;
#endif

        execute_command(command);

#ifdef SUANPAN_WIN
        SetConsoleTextAttribute(handle, current_attribute);
#else
        SUANPAN_SYNC_COUT << FOREGROUND_GREEN;
#endif

        return SUANPAN_SUCCESS;
    }

    if(is_equal(command_id, "version")) print_version();
    else suanpan_error("command not found.\n");

    return SUANPAN_SUCCESS;
}

int process_file(const shared_ptr<Bead>& model, const char* file_name) {
    ifstream input_file;

    input_file.open(fs::path(file_name));

    if(!input_file.is_open()) {
        string new_name = file_name;
        new_name += ".supan";
        input_file.open(fs::path(new_name));
        if(!input_file.is_open()) {
            new_name = file_name;
            suanpan::to_upper(new_name);
            new_name += ".supan";
            input_file.open(fs::path(new_name));
            if(!input_file.is_open()) {
                new_name = file_name;
                suanpan::to_lower(new_name);
                new_name += ".supan";
                input_file.open(fs::path(new_name));
                if(!input_file.is_open()) {
                    suanpan_error("process_file() cannot open %s.\n", fs::path(file_name).generic_string().c_str());
                    return SUANPAN_EXIT;
                }
            }
        }
    }

    string all_line, command_line;
    while(!getline(input_file, command_line).fail())
        if(!command_line.empty() && command_line[0] != '#' && command_line[0] != '!') {
            if(const auto if_comment = command_line.find('!'); string::npos != if_comment) command_line.erase(if_comment);
            for(auto& c : command_line) if(',' == c || '\t' == c || '\r' == c || '\n' == c) c = ' ';
            while(*command_line.crbegin() == ' ') command_line.pop_back();
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

int process_file(const shared_ptr<Bead>& model, istringstream& command) {
    string file_name;
    if(!get_input(command, file_name)) {
        suanpan_error("process_file() needs a file name.\n");
        return SUANPAN_SUCCESS;
    }

    return process_file(model, file_name.c_str());
}

int create_new_domain(const shared_ptr<Bead>& model, istringstream& command) {
    unsigned domain_id;
    if(!get_input(command, domain_id)) {
        suanpan_error("create_new_domain() requires a tag.\n");
        return SUANPAN_SUCCESS;
    }

    model->set_current_domain_tag(domain_id);

    if(auto& tmp_domain = get_domain(model, domain_id); nullptr != tmp_domain) suanpan_info("create_new_domain() switches to Domain %u.\n", domain_id);
    else {
        tmp_domain = make_shared<Domain>(domain_id);
        if(nullptr != tmp_domain) suanpan_info("create_new_domain() successfully creates Domain %u.\n", domain_id);
    }

    return SUANPAN_SUCCESS;
}

int disable_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(nullptr == domain) {
        suanpan_error("disable_object() needs a valid domain.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("disable_object() needs object type.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "domain")) while(get_input(command, tag)) model->disable_domain(tag);
    else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->disable_step(tag);
    else if(is_equal(object_type, "converger")) while(get_input(command, tag)) domain->disable_converger(tag);
    else if(is_equal(object_type, "constraint")) while(get_input(command, tag)) domain->disable_constraint(tag);
    else if(is_equal(object_type, "element")) while(get_input(command, tag)) domain->disable_element(tag);
    else if(is_equal(object_type, "load")) while(get_input(command, tag)) domain->disable_load(tag);
    else if(is_equal(object_type, "material")) while(get_input(command, tag)) domain->disable_material(tag);
    else if(is_equal(object_type, "node")) while(get_input(command, tag)) domain->disable_node(tag);
    else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->disable_recorder(tag);
    else if(is_equal(object_type, "modifier")) while(get_input(command, tag)) domain->disable_modifier(tag);
    else if(is_equal(object_type, "print")) SUANPAN_PRINT = false;

    return SUANPAN_SUCCESS;
}

int enable_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(nullptr == domain) {
        suanpan_error("enable_object() needs a valid domain.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("enable_object() needs object type.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "domain")) while(get_input(command, tag)) model->enable_domain(tag);
    else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->enable_step(tag);
    else if(is_equal(object_type, "converger")) while(get_input(command, tag)) domain->enable_converger(tag);
    else if(is_equal(object_type, "constraint")) while(get_input(command, tag)) domain->enable_constraint(tag);
    else if(is_equal(object_type, "element")) while(get_input(command, tag)) domain->enable_element(tag);
    else if(is_equal(object_type, "load")) while(get_input(command, tag)) domain->enable_load(tag);
    else if(is_equal(object_type, "material")) while(get_input(command, tag)) domain->enable_material(tag);
    else if(is_equal(object_type, "node")) while(get_input(command, tag)) domain->enable_node(tag);
    else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->enable_recorder(tag);
    else if(is_equal(object_type, "modifier")) while(get_input(command, tag)) domain->enable_modifier(tag);
    else if(is_equal(object_type, "print")) SUANPAN_PRINT = true;

    return SUANPAN_SUCCESS;
}

int erase_object(const shared_ptr<Bead>& model, istringstream& command) {
    const auto& domain = get_current_domain(model);
    if(nullptr == domain) {
        suanpan_error("erase_object() needs a valid domain.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("erase_object() needs object type.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "domain")) while(get_input(command, tag)) model->erase_domain(tag);
    else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->erase_step(tag);
    else if(is_equal(object_type, "converger")) while(get_input(command, tag)) domain->erase_converger(tag);
    else if(is_equal(object_type, "constraint")) while(get_input(command, tag)) domain->erase_constraint(tag);
    else if(is_equal(object_type, "element")) while(get_input(command, tag)) domain->erase_element(tag);
    else if(is_equal(object_type, "load")) while(get_input(command, tag)) domain->erase_load(tag);
    else if(is_equal(object_type, "material")) while(get_input(command, tag)) domain->erase_material(tag);
    else if(is_equal(object_type, "node")) while(get_input(command, tag)) domain->erase_node(tag);
    else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->erase_recorder(tag);
    else if(is_equal(object_type, "modifier")) while(get_input(command, tag)) domain->erase_modifier(tag);

    return SUANPAN_SUCCESS;
}

int save_object(const shared_ptr<DomainBase>& domain, istringstream& command) {
    if(nullptr == domain) {
        suanpan_error("erase_object() needs a valid domain.\n");
        return SUANPAN_SUCCESS;
    }

    string object_id;
    if(!get_input(command, object_id)) {
        suanpan_error("save_object() needs a valid object type.\n");
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
        suanpan_error("list_object() needs a valid domain.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("list_object() needs object type.\n");
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

    suanpan_info("This domain has the following %ss:", object_type.c_str());
    for(const auto& I : list) suanpan_info("\t%u", I);
    suanpan_info(".\n");

    return SUANPAN_SUCCESS;
}

int suspend_object(const shared_ptr<DomainBase>& domain, istringstream& command) {
    if(nullptr == domain) {
        suanpan_error("suspend_object() needs a valid domain.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("suspend_object() needs object type.\n");
        return SUANPAN_SUCCESS;
    }

    const auto step_tag = domain->get_current_step_tag();

    if(unsigned tag; is_equal(object_type, "constraint")) while(!command.eof() && get_input(command, tag)) { if(domain->find_constraint(tag)) domain->get_constraint(tag)->set_end_step(step_tag); }
    else if(is_equal(object_type, "load")) while(!command.eof() && get_input(command, tag)) { if(domain->find_load(tag)) domain->get_load(tag)->set_end_step(step_tag); }

    return SUANPAN_SUCCESS;
}

int protect_object(const shared_ptr<DomainBase>& domain, istringstream& command) {
    if(nullptr == domain) {
        suanpan_error("protect_object() needs a valid domain.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("protect_object() needs object type.\n");
        return SUANPAN_SUCCESS;
    }

    if(unsigned tag; is_equal(object_type, "element")) while(!command.eof() && get_input(command, tag)) { if(domain->find<Element>(tag)) domain->get<Element>(tag)->guard(); }
    else if(is_equal(object_type, "node")) while(!command.eof() && get_input(command, tag)) { if(domain->find<Node>(tag)) domain->get<Node>(tag)->guard(); }

    return SUANPAN_SUCCESS;
}

int create_new_acceleration(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_acceleration() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_acceleration() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_acceleration() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_acceleration() needs dof.\n");
        return SUANPAN_SUCCESS;
    }

    uword node_id;
    vector<uword> node_pool;
    while(get_input(command, node_id)) node_pool.emplace_back(node_id);

    if(const auto step_tag = domain->get_current_step_tag(); !domain->insert(make_shared<NodalAcceleration>(load_id, step_tag, magnitude, uvec(node_pool), dof_id, amplitude_id))) suanpan_error("create_new_acceleration() fails to create new load.\n");

    return SUANPAN_SUCCESS;
}

int create_new_amplitude(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string amplitude_type;
    if(!get_input(command, amplitude_type)) {
        suanpan_error("create_new_amplitude() needs a valid amplitude type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_amplitude() needs a valid amplitude type.\n");
        return SUANPAN_SUCCESS;
    }

    if(const auto step_tag = domain->get_current_step_tag(); is_equal(amplitude_type, "Constant")) domain->insert(make_shared<Constant>(tag, step_tag));
    else if(is_equal(amplitude_type, "Ramp")) domain->insert(make_shared<Ramp>(tag, step_tag));
    else if(is_equal(amplitude_type, "Tabular")) {
        string file_name;
        if(!get_input(command, file_name)) {
            suanpan_error("create_new_amplitude() needs a valid file.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(make_shared<Tabular>(tag, std::move(file_name), step_tag));
    }
    else if(is_equal(amplitude_type, "Decay")) {
        double A;
        if(!get_input(command, A)) {
            suanpan_error("create_new_amplitude() needs a A.\n");
            return SUANPAN_SUCCESS;
        }
        double TD;
        if(!get_input(command, TD)) {
            suanpan_error("create_new_amplitude() needs a TD.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(make_shared<Decay>(tag, A, TD, step_tag));
    }
    else if(is_equal(amplitude_type, "Linear")) {
        double A;
        if(!get_input(command, A)) {
            suanpan_error("create_new_amplitude() needs a slope.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(make_shared<Linear>(tag, A, step_tag));
    }
    else if(is_equal(amplitude_type, "Combine")) {
        vector<uword> tag_pool;
        uword t_tag;
        while(get_input(command, t_tag)) tag_pool.emplace_back(t_tag);
        domain->insert(make_shared<Combine>(tag, uvec(tag_pool), step_tag));
    }
    else if(is_equal(amplitude_type, "Modulated") || is_equal(amplitude_type, "Sine") || is_equal(amplitude_type, "Cosine")) {
        double W;
        if(!get_input(command, W)) {
            suanpan_error("create_new_amplitude() needs a period/amplitude.\n");
            return SUANPAN_SUCCESS;
        }

        double amp;
        vector<double> A;
        while(get_input(command, amp)) A.emplace_back(amp);

        if(is_equal(amplitude_type, "Modulated")) domain->insert(make_shared<Modulated>(tag, W, std::move(A), step_tag));
        else if(is_equal(amplitude_type, "Sine")) domain->insert(make_shared<Sine>(tag, W, std::move(A), step_tag));
        else if(is_equal(amplitude_type, "Cosine")) domain->insert(make_shared<Cosine>(tag, W, std::move(A), step_tag));
    }
    else if(is_equal(amplitude_type, "NZStrongMotion")) {
        string name;
        if(!get_input(command, name)) {
            suanpan_error("create_new_amplitude() needs a name.\n");
            return SUANPAN_SUCCESS;
        }

        domain->insert(make_shared<NZStrongMotion>(tag, name.c_str(), step_tag));
    }

    return SUANPAN_SUCCESS;
}

int create_new_bc(const shared_ptr<DomainBase>& domain, istringstream& command, const bool flag) {
    unsigned bc_id;
    if(!get_input(command, bc_id)) {
        suanpan_error("create_new_bc() needs bc tag.\n");
        return SUANPAN_SUCCESS;
    }

    string dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_bc() needs valid DoFs.\n");
        return SUANPAN_SUCCESS;
    }

    uword node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    const auto bc_type = suanpan::to_lower(dof_id[0]);

    if(const auto step_tag = domain->get_current_step_tag(); flag) {
        if(is_equal(bc_type, 'p')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), "PINNED"));
        else if(is_equal(bc_type, 'e')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), "ENCASTRE"));
        else if(is_equal(bc_type, 'x')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), "XSYMM"));
        else if(is_equal(bc_type, 'y')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), "YSYMM"));
        else if(is_equal(bc_type, 'z')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), "ZSYMM"));
        else if(is_equal(bc_type, '1')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), 1));
        else if(is_equal(bc_type, '2')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), 2));
        else if(is_equal(bc_type, '3')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), 3));
        else if(is_equal(bc_type, '4')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), 4));
        else if(is_equal(bc_type, '5')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), 5));
        else if(is_equal(bc_type, '6')) domain->insert(make_shared<PenaltyBC>(bc_id, step_tag, uvec(node_tag), 6));
    }
    else {
        if(is_equal(bc_type, 'p')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), "PINNED"));
        else if(is_equal(bc_type, 'e')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), "ENCASTRE"));
        else if(is_equal(bc_type, 'x')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), "XSYMM"));
        else if(is_equal(bc_type, 'y')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), "YSYMM"));
        else if(is_equal(bc_type, 'z')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), "ZSYMM"));
        else if(is_equal(bc_type, '1')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), 1));
        else if(is_equal(bc_type, '2')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), 2));
        else if(is_equal(bc_type, '3')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), 3));
        else if(is_equal(bc_type, '4')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), 4));
        else if(is_equal(bc_type, '5')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), 5));
        else if(is_equal(bc_type, '6')) domain->insert(make_shared<MultiplierBC>(bc_id, step_tag, uvec(node_tag), 6));
    }

    return SUANPAN_SUCCESS;
}

int create_new_groupbc(const shared_ptr<DomainBase>& domain, istringstream& command, const bool flag) {
    unsigned bc_id;
    if(!get_input(command, bc_id)) {
        suanpan_error("create_new_groupbc() needs bc tag.\n");
        return SUANPAN_SUCCESS;
    }

    string dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_groupbc() needs valid DoFs.\n");
        return SUANPAN_SUCCESS;
    }

    uword group;
    vector<uword> group_tag;
    while(get_input(command, group)) group_tag.push_back(group);

    const auto bc_type = suanpan::to_lower(dof_id[0]);

    if(const auto step_tag = domain->get_current_step_tag(); flag) {
        if(is_equal(bc_type, 'p')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), "PINNED"));
        else if(is_equal(bc_type, 'e')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), "ENCASTRE"));
        else if(is_equal(bc_type, 'x')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), "XSYMM"));
        else if(is_equal(bc_type, 'y')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), "YSYMM"));
        else if(is_equal(bc_type, 'z')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), "ZSYMM"));
        else if(is_equal(bc_type, '1')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), 1));
        else if(is_equal(bc_type, '2')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), 2));
        else if(is_equal(bc_type, '3')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), 3));
        else if(is_equal(bc_type, '4')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), 4));
        else if(is_equal(bc_type, '5')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), 5));
        else if(is_equal(bc_type, '6')) domain->insert(make_shared<GroupPenaltyBC>(bc_id, step_tag, uvec(group_tag), 6));
    }
    else {
        if(is_equal(bc_type, 'p')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), "PINNED"));
        else if(is_equal(bc_type, 'e')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), "ENCASTRE"));
        else if(is_equal(bc_type, 'x')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), "XSYMM"));
        else if(is_equal(bc_type, 'y')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), "YSYMM"));
        else if(is_equal(bc_type, 'z')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), "ZSYMM"));
        else if(is_equal(bc_type, '1')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), 1));
        else if(is_equal(bc_type, '2')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), 2));
        else if(is_equal(bc_type, '3')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), 3));
        else if(is_equal(bc_type, '4')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), 4));
        else if(is_equal(bc_type, '5')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), 5));
        else if(is_equal(bc_type, '6')) domain->insert(make_shared<GroupMultiplierBC>(bc_id, step_tag, uvec(group_tag), 6));
    }

    return SUANPAN_SUCCESS;
}

int create_new_bodyforce(const shared_ptr<DomainBase>& domain, istringstream& command, const bool flag) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_bodyforce() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_bodyforce() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_bodyforce() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_bodyforce() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned element;
    vector<uword> element_tag;
    while(get_input(command, element)) element_tag.push_back(element);

    if(flag) { if(!domain->insert(make_shared<GroupBodyForce>(load_id, domain->get_current_step_tag(), magnitude, uvec(element_tag), dof_id, amplitude_id))) suanpan_error("create_new_bodyforce() fails to create new load.\n"); }
    else if(!domain->insert(make_shared<BodyForce>(load_id, domain->get_current_step_tag(), magnitude, uvec(element_tag), dof_id, amplitude_id))) suanpan_error("create_new_bodyforce() fails to create new load.\n");

    return SUANPAN_SUCCESS;
}

int create_new_cload(const shared_ptr<DomainBase>& domain, istringstream& command, const bool flag) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_cload() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_cload() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_cload() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_cload() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(flag) { if(!domain->insert(make_shared<GroupNodalForce>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_cload() fails to create new load.\n"); }
    else { if(!domain->insert(make_shared<NodalForce>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_cload() fails to create new load.\n"); }

    return SUANPAN_SUCCESS;
}

int create_new_lineudl(const shared_ptr<DomainBase>& domain, istringstream& command, const unsigned dimension) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_lineudl() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_lineudl() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_lineudl() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_lineudl() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(2 == dimension) { if(!domain->insert(make_shared<LineUDL2D>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_lineudl() fails to create new load.\n"); }
    else { if(!domain->insert(make_shared<LineUDL3D>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_lineudl() fails to create new load.\n"); }

    return SUANPAN_SUCCESS;
}

int create_new_converger(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string converger_id;
    if(!get_input(command, converger_id)) {
        suanpan_error("create_new_converger() requires converger type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_converger() requires a tag.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    if(is_equal(converger_id.substr(0, 5), "Logic")) {
        unsigned tag_a, tag_b;
        if(!get_input(command, tag_a) || !get_input(command, tag_b)) {
            suanpan_error("create_new_converger() requires a tag.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(converger_id, "LogicAND") && domain->insert(make_shared<LogicAND>(tag, tag_a, tag_b))) code = 1; // NOLINT(bugprone-branch-clone)
        else if(is_equal(converger_id, "LogicOR") && domain->insert(make_shared<LogicOR>(tag, tag_a, tag_b))) code = 1;
        else if(is_equal(converger_id, "LogicXOR") && domain->insert(make_shared<LogicXOR>(tag, tag_a, tag_b))) code = 1;
        else suanpan_error("create_new_converger() cannot identify the converger type.\n");
    }
    else {
        auto tolerance = 1E-6;
        if(!is_equal(converger_id, "FixedNumber") && (!command.eof() && !get_input(command, tolerance))) {
            suanpan_error("create_new_converger() reads wrong tolerance.\n");
            return SUANPAN_SUCCESS;
        }

        auto max_iteration = 10;
        if(!command.eof() && !get_input(command, max_iteration)) {
            suanpan_error("create_new_converger() reads wrong max iteration.\n");
            return SUANPAN_SUCCESS;
        }

        string print_flag = "false";
        if(!command.eof() && !get_input(command, print_flag)) {
            suanpan_error("create_new_converger() reads wrong print flag.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(converger_id, "AbsResidual") && domain->insert(make_shared<AbsResidual>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1; // NOLINT(bugprone-branch-clone)
        else if(is_equal(converger_id, "RelResidual") && domain->insert(make_shared<RelResidual>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsIncreDisp") && domain->insert(make_shared<AbsIncreDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelIncreDisp") && domain->insert(make_shared<RelIncreDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsDisp") && domain->insert(make_shared<AbsDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelDisp") && domain->insert(make_shared<RelDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsError") && domain->insert(make_shared<AbsError>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelError") && domain->insert(make_shared<RelError>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsIncreEnergy") && domain->insert(make_shared<AbsIncreEnergy>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelIncreEnergy") && domain->insert(make_shared<RelIncreEnergy>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "FixedNumber") && domain->insert(make_shared<FixedNumber>(tag, max_iteration, is_true(print_flag)))) code = 1;
        else suanpan_error("create_new_converger() cannot identify the converger type.\n");
    }

    if(1 == code) {
        if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_converger_tag(tag);
        domain->set_current_converger_tag(tag);
    }
    else suanpan_error("create_new_converger() fails to create the new converger.\n");

    return SUANPAN_SUCCESS;
}

int create_new_criterion(const shared_ptr<DomainBase>& domain, istringstream& command) {
    const auto& step_tag = domain->get_current_step_tag();
    if(0 == step_tag) {
        suanpan_error("create_new_criterion() needs a valid step.\n");
        return SUANPAN_SUCCESS;
    }

    string criterion_type;
    if(!get_input(command, criterion_type)) {
        suanpan_error("create_new_criterion() need a criterion type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_criterion() requires a tag.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(criterion_type.substr(0, 5), "Logic")) {
        unsigned tag_a, tag_b;
        if(!get_input(command, tag_a) || !get_input(command, tag_b)) {
            suanpan_error("create_new_criterion() requires a valid tag.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(criterion_type, "LogicCriterionAND")) domain->insert(make_shared<LogicCriterionAND>(tag, step_tag, tag_a, tag_b));
        else if(is_equal(criterion_type, "LogicCriterionOR")) domain->insert(make_shared<LogicCriterionOR>(tag, step_tag, tag_a, tag_b));

        return SUANPAN_SUCCESS;
    }

    if(is_equal(criterion_type, "StrainEnergyEvolution")) {
        unsigned incre_level, final_level;
        if(!get_input(command, incre_level)) {
            suanpan_error("create_new_criterion() requires a valid level.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, final_level)) {
            suanpan_error("create_new_criterion() requires a valid level.\n");
            return SUANPAN_SUCCESS;
        }

        auto weight = 1.;
        if(!command.eof() && !get_input(command, weight)) {
            suanpan_error("create_new_criterion() requires a valid weight of central element.\n");
            return SUANPAN_SUCCESS;
        }
        auto iteration = 2;
        if(!command.eof() && !get_input(command, iteration)) {
            suanpan_error("create_new_criterion() requires a valid number of iteration.\n");
            return SUANPAN_SUCCESS;
        }
        auto reactivation = 10;
        if(!command.eof() && !get_input(command, reactivation)) {
            suanpan_error("create_new_criterion() requires a valid number of reactivation ratio.\n");
            return SUANPAN_SUCCESS;
        }
        auto propagation = .5;
        if(!command.eof() && !get_input(command, propagation)) {
            suanpan_error("create_new_criterion() requires a valid propagation factor.\n");
            return SUANPAN_SUCCESS;
        }
        auto tolerance = 1E-5;
        if(!command.eof() && !get_input(command, tolerance)) {
            suanpan_error("create_new_criterion() requires a valid tolerance.\n");
            return SUANPAN_SUCCESS;
        }

        domain->insert(make_shared<StrainEnergyEvolution>(tag, step_tag, incre_level, final_level, weight, iteration, reactivation, propagation, tolerance));

        return SUANPAN_SUCCESS;
    }

    if(is_equal(criterion_type, "MaxHistory")) {
        string type;
        double limit;
        if(!get_input(command, type)) {
            suanpan_error("create_new_criterion() requires a valid type.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, limit)) {
            suanpan_error("create_new_criterion() requires a valid limit.\n");
            return SUANPAN_SUCCESS;
        }

        domain->insert(make_shared<MaxHistory>(tag, step_tag, to_list(type.c_str()), limit));

        return SUANPAN_SUCCESS;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("create_new_criterion() requires a node.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof;
    if(!get_input(command, dof)) {
        suanpan_error("create_new_criterion() requires a dof.\n");
        return SUANPAN_SUCCESS;
    }

    double limit;
    if(!get_input(command, limit)) {
        suanpan_error("create_new_criterion() requires a limit.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(criterion_type, "MaxDisplacement")) domain->insert(make_shared<MaxDisplacement>(tag, step_tag, node, dof, limit));
    else if(is_equal(criterion_type, "MinDisplacement")) domain->insert(make_shared<MinDisplacement>(tag, step_tag, node, dof, limit));
    else if(is_equal(criterion_type, "MaxResistance")) domain->insert(make_shared<MaxResistance>(tag, step_tag, node, dof, limit));
    else if(is_equal(criterion_type, "MinResistance")) domain->insert(make_shared<MinResistance>(tag, step_tag, node, dof, limit));

    return SUANPAN_SUCCESS;
}

int create_new_displacement(const shared_ptr<DomainBase>& domain, istringstream& command, const bool flag) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_displacement() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_displacement() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_displacement() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_displacement() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(const auto step_tag = domain->get_current_step_tag(); flag) { if(!domain->insert(make_shared<GroupNodalDisplacement>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_displacement() fails to create new load.\n"); }
    else { if(!domain->insert(make_shared<NodalDisplacement>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_displacement() fails to create new load.\n"); }

    return SUANPAN_SUCCESS;
}

int create_new_nodegroup(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_nodegroup() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    uword value;
    vector<uword> value_pool;
    while(get_input(command, value)) value_pool.push_back(value);

    if(!domain->insert(make_shared<NodeGroup>(tag, value_pool))) suanpan_error("create_new_nodegroup() fails to create new node group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_elementgroup(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_elementgroup() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    uword value;
    vector<uword> value_pool;
    while(get_input(command, value)) value_pool.push_back(value);

    if(!domain->insert(make_shared<ElementGroup>(tag, value_pool))) suanpan_error("create_new_elementgroup() fails to create new element group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_generate(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("create_new_generate() needs a type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_generate() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    int start, interval, end;
    if(!get_input(command, start)) {
        suanpan_error("create_new_generate() needs a valid tag.\n");
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

    if(is_equal(type, "nodegroup") && !domain->insert(make_shared<NodeGroup>(tag, tag_pool))) suanpan_error("create_new_generate() fails to create new node group.\n");
    else if(is_equal(type, "elementgroup") && !domain->insert(make_shared<ElementGroup>(tag, tag_pool))) suanpan_error("create_new_generate() fails to create new element group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_generatebyrule(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("create_new_generatebyrule() needs a type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_generatebyrule() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof;
    if(!get_input(command, dof)) {
        suanpan_error("create_new_generatebyrule() needs a valid dof.\n");
        return SUANPAN_SUCCESS;
    }

    double para;
    vector<double> pool;
    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(is_equal(type, "nodegroup") && !domain->insert(make_shared<NodeGroup>(tag, dof, pool))) suanpan_error("create_new_generatebyrule() fails to create new node group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_generatebyplane(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("create_new_generatebyplane() needs a type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_generatebyplane() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    double para;
    vector<double> pool;
    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(pool.empty()) return SUANPAN_SUCCESS;

    if(is_equal(type, "nodegroup") && !domain->insert(make_shared<NodeGroup>(tag, pool))) suanpan_error("create_new_generatebyplane() fails to create new node group.\n");

    return SUANPAN_SUCCESS;
}

int create_new_generatebypoint(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("create_new_generatebypoint() needs a type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_generatebypoint() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    double para;
    vector<double> pool;
    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(pool.size() % 2 == 0) { if(const auto size = static_cast<long long>(pool.size()) / 2; is_equal(type, "nodegroup") && !domain->insert(make_shared<NodeGroup>(tag, vector(pool.begin(), pool.begin() + size), vector(pool.end() - size, pool.end())))) suanpan_error("create_new_generatebypoint() fails to create new node group.\n"); }

    return SUANPAN_SUCCESS;
}

int create_new_groupgroup(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_groupgroup() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    uword para;
    vector<uword> pool;
    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(!domain->insert(make_shared<GroupGroup>(tag, pool))) suanpan_error("create_new_groupgroup() fails to create new group of groups.\n");

    return SUANPAN_SUCCESS;
}

int create_new_external_module(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string library_name;

    if(!get_input(command, library_name)) {
        suanpan_error("create_new_external_module() needs module name.\n");
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

int create_new_fixedlength(const shared_ptr<DomainBase>& domain, istringstream& command, const unsigned dof) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_fixedlength() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    uword node_i, node_j;
    if(!get_input(command, node_i) || !get_input(command, node_j)) {
        suanpan_error("create_new_fixedlength() needs two node tags.\n");
        return SUANPAN_SUCCESS;
    }

    domain->insert(make_unique<FixedLength>(tag, domain->get_current_step_tag(), dof, uvec{node_i, node_j}));

    return SUANPAN_SUCCESS;
}

int create_new_initial(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string variable_type;
    if(!get_input(command, variable_type)) {
        suanpan_error("create_new_initial() needs a valid variable type.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal("material", variable_type)) {
        string state_type;
        if(!get_input(command, state_type)) {
            suanpan_error("create_new_initial() needs a valid state type.\n");
            return SUANPAN_SUCCESS;
        }

        unsigned mat_tag;
        if(!get_input(command, mat_tag)) {
            suanpan_error("create_new_initial() needs a valid material tag.\n");
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
                suanpan_error("create_new_initial() needs a valid magnitude.\n");
                return SUANPAN_SUCCESS;
            }

        unsigned ref_node;
        if(!get_input(command, ref_node) || !domain->find_node(ref_node)) {
            suanpan_error("create_new_initial() needs a valid reference node tag.\n");
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
                suanpan_error("create_new_initial() needs a valid node tag.\n");
                return SUANPAN_SUCCESS;
            }
        }

        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_initial() needs a valid magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_tag;
    if(!get_input(command, dof_tag)) {
        suanpan_error("create_new_initial() needs a valid dof tag.\n");
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
                suanpan_error("create_new_initial() needs a valid node tag.\n");
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
                suanpan_error("create_new_initial() needs a valid node tag.\n");
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
                suanpan_error("create_new_initial() needs a valid node tag.\n");
                return SUANPAN_SUCCESS;
            }
        }

    return SUANPAN_SUCCESS;
}

int create_new_integrator(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string integrator_type;
    if(!get_input(command, integrator_type)) {
        suanpan_error("create_new_integrator() needs a valid integrator type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_integrator() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    if(if_contain(suanpan::to_upper(cref(integrator_type)), suanpan::to_upper(string("Newmark")))) {
        auto alpha = .25, beta = .5;
        if(!command.eof()) {
            if(!get_input(command, alpha)) {
                suanpan_error("create_new_integrator() needs a valid alpha.\n");
                return SUANPAN_SUCCESS;
            }
            if(!get_input(command, beta)) {
                suanpan_error("create_new_integrator() needs a valid beta.\n");
                return SUANPAN_SUCCESS;
            }
        }

        if(is_equal(integrator_type, "Newmark")) { if(domain->insert(make_shared<Newmark>(tag, alpha, beta))) code = 1; }
        else if(is_equal(integrator_type, "RayleighNewmark")) {
            vec p(4, fill::zeros);
            auto idx = 0llu;
            while(!command.eof() && idx < p.n_elem)
                if(!get_input(command, p(idx++))) {
                    suanpan_error("create_new_integrator() needs a valid parameter for Rayleigh damping.\n");
                    return SUANPAN_SUCCESS;
                }

            if(domain->insert(make_shared<RayleighNewmark>(tag, alpha, beta, p(0), p(1), p(2), p(3)))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeNewmark")) {
            vector<double> damping_coef, frequency;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("create_new_integrator() needs a valid damping coefficient.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("create_new_integrator() needs a valid frequency.\n");
                    return SUANPAN_SUCCESS;
                }
                frequency.emplace_back(t_para);
            }

            if(domain->insert(make_shared<LeeNewmark>(tag, damping_coef, frequency, alpha, beta))) code = 1;
        }
        else if(integrator_type.size() >= 14 && is_equal(integrator_type.substr(0, 14), "LeeNewmarkFull")) {
            vector<LeeNewmarkFull::Mode> modes;

            auto omega = 0., zeta = 0., para_a = .0, para_b = .0;

            auto get_basic_input = [&] {
                if(!get_input(command, zeta)) {
                    suanpan_error("create_new_integrator() needs a valid zeta_p.\n");
                    return SUANPAN_FAIL;
                }
                if(!get_input(command, omega)) {
                    suanpan_error("create_new_integrator() needs a valid omega_p.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            auto get_first = [&] {
                if(!get_input(command, para_a)) {
                    suanpan_error("create_new_integrator() needs a valid parameter.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            auto get_second = [&] {
                if(!get_input(command, para_b)) {
                    suanpan_error("create_new_integrator() needs a valid parameter.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            while(!command.eof()) {
                string type;
                if(!get_input(command, type)) {
                    suanpan_error("create_new_integrator() needs a valid type.\n");
                    return SUANPAN_SUCCESS;
                }
                if(is_equal("-type0", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T0, vec{}, zeta, omega});
                }
                else if(is_equal("-type1", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T1, vec{static_cast<double>(static_cast<unsigned>(para_a))}, zeta, omega});
                }
                else if(is_equal("-type2", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T2, vec{static_cast<double>(static_cast<unsigned>(para_a)), static_cast<double>(static_cast<unsigned>(para_b))}, zeta, omega});
                }
                else if(is_equal("-type3", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T3, vec{para_a}, zeta, omega});
                }
                else if(is_equal("-type4", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return SUANPAN_SUCCESS;
                    double para_c, para_d, para_e;
                    if(!get_input(command, para_c) || !get_input(command, para_d) || !get_input(command, para_e)) {
                        suanpan_error("create_new_integrator() needs a valid parameter.\n");
                        return SUANPAN_SUCCESS;
                    }
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T4, vec{static_cast<double>(static_cast<unsigned>(para_a)), static_cast<double>(static_cast<unsigned>(para_b)), static_cast<double>(static_cast<unsigned>(para_c)), static_cast<double>(static_cast<unsigned>(para_d)), para_e}, zeta, omega});
                }
                else {
                    suanpan_error("create_new_integrator() needs a valid type.\n");
                    return SUANPAN_SUCCESS;
                }
            }

            if(is_equal(integrator_type, "LeeNewmarkFullTrial")) { if(domain->insert(make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::TRIAL))) code = 1; }
            else if(is_equal(integrator_type, "LeeNewmarkFullCurrent") || is_equal(integrator_type, "LeeNewmarkFull")) { if(domain->insert(make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::CURRENT))) code = 1; }
            else if(is_equal(integrator_type, "LeeNewmarkFullInitial")) { if(domain->insert(make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::INITIAL))) code = 1; }
        }
        else if(is_equal(integrator_type, "WilsonPenzienNewmark")) {
            vector<double> damping_coef;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("create_new_integrator() needs a valid damping coefficient.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
            }

            if(domain->insert(make_shared<WilsonPenzienNewmark>(tag, damping_coef, alpha, beta))) code = 1;
        }
    }
    else if(is_equal(integrator_type, "GeneralizedAlpha") || is_equal(integrator_type, "GeneralisedAlpha")) {
        vector<double> pool;
        pool.reserve(2);

        double para;
        while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

        if(pool.empty() && domain->insert(make_shared<GeneralizedAlpha>(tag, .5))) code = 1; // NOLINT(bugprone-branch-clone)
        else if(1 == pool.size() && domain->insert(make_shared<GeneralizedAlpha>(tag, std::min(std::max(0., pool[0]), 1.)))) code = 1;
        else if(2 == pool.size() && domain->insert(make_shared<GeneralizedAlpha>(tag, pool[0], pool[1]))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSU0")) {
        vec pool(3);

        for(auto& I : pool)
            if(!get_input(command, I)) {
                suanpan_error("create_new_integrator() needs a valid damping radius.\n");
                return SUANPAN_SUCCESS;
            }

        if(domain->insert(make_shared<GSSSSU0>(tag, std::move(pool)))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSV0")) {
        vec pool(3);

        for(auto& I : pool)
            if(!get_input(command, I)) {
                suanpan_error("create_new_integrator() needs a valid damping radius.\n");
                return SUANPAN_SUCCESS;
            }

        if(domain->insert(make_shared<GSSSSV0>(tag, std::move(pool)))) code = 1;
    }
    else if(is_equal(integrator_type, "BatheTwoStep") && domain->insert(make_shared<BatheTwoStep>(tag))) code = 1;

    if(1 == code) {
        if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_integrator_tag(tag);
        domain->set_current_integrator_tag(tag);
    }
    else suanpan_error("create_new_integrator() fails to create the new integrator.\n");

    return SUANPAN_SUCCESS;
}

int create_new_mass(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_mass() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("create_new_mass() needs one valid node.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_mass() needs a valid magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof;
    vector<uword> dof_tag;
    while(get_input(command, dof)) dof_tag.push_back(dof);

    domain->insert(make_shared<Mass>(tag, node, magnitude, uvec(dof_tag)));

    return SUANPAN_SUCCESS;
}

int create_new_modifier(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string modifier_type;
    if(!get_input(command, modifier_type)) {
        suanpan_error("create_new_modifier() needs a valid modifier type.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Modifier> new_modifier = nullptr;

    if(is_equal(modifier_type, "LumpedSimple")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("create_new_modifier() needs a valid tag.\n");
            return SUANPAN_SUCCESS;
        }

        vector<uword> element_tag;
        unsigned e_tag;
        while(!command.eof() && get_input(command, e_tag)) element_tag.emplace_back(e_tag);

        new_modifier = make_unique<LumpedSimple>(tag, uvec(element_tag));
    }
    else if(is_equal(modifier_type, "LumpedScale")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("create_new_modifier() needs a valid tag.\n");
            return SUANPAN_SUCCESS;
        }

        vector<uword> element_tag;
        unsigned e_tag;
        while(!command.eof() && get_input(command, e_tag)) element_tag.emplace_back(e_tag);

        new_modifier = make_unique<LumpedScale>(tag, uvec(element_tag));
    }
    else if(is_equal(modifier_type, "Rayleigh")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("create_new_modifier() needs a valid tag.\n");
            return SUANPAN_SUCCESS;
        }

        double a, b, c, d;
        if(!get_input(command, a)) {
            suanpan_error("create_new_modifier() needs four valid numbers.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, b)) {
            suanpan_error("create_new_modifier() needs four valid numbers.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, c)) {
            suanpan_error("create_new_modifier() needs four valid numbers.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, d)) {
            suanpan_error("create_new_modifier() needs four valid numbers.\n");
            return SUANPAN_SUCCESS;
        }
        vector<uword> element_tag;
        unsigned e_tag;
        while(!command.eof() && get_input(command, e_tag)) element_tag.emplace_back(e_tag);

        new_modifier = make_unique<Rayleigh>(tag, a, b, c, d, uvec(element_tag));
    }
    else if(is_equal(modifier_type, "ElementalModal")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("create_new_modifier() needs a valid tag.\n");
            return SUANPAN_SUCCESS;
        }

        double a, b;
        if(!get_input(command, a)) {
            suanpan_error("create_new_modifier() needs two valid numbers.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, b)) {
            suanpan_error("create_new_modifier() needs two valid numbers.\n");
            return SUANPAN_SUCCESS;
        }

        vector<uword> element_tag;
        unsigned e_tag;
        while(!command.eof() && get_input(command, e_tag)) element_tag.emplace_back(e_tag);

        new_modifier = make_unique<ElementalModal>(tag, a, b, uvec(element_tag));
    }
    else {
        // check if the library is already loaded
        auto code = false;
        for(const auto& I : domain->get_external_module_pool())
            if(is_equal(I->library_name, modifier_type) || I->locate_cpp_module(modifier_type)) {
                code = true;
                break;
            }

        // not loaded then try load it
        if(!code && domain->insert(make_shared<ExternalModule>(modifier_type))) code = true;

        // if loaded find corresponding function
        if(code)
            for(const auto& I : domain->get_external_module_pool()) {
                if(I->locate_cpp_module(modifier_type)) I->new_object(new_modifier, command);
                if(new_modifier != nullptr) break;
            }
    }

    if(nullptr == new_modifier || !domain->insert(std::move(new_modifier))) suanpan_error("create_new_modifier() fails to create new modifier.\n");

    return SUANPAN_SUCCESS;
}

int create_new_mpc(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_mpc() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude;
    if(!get_input(command, amplitude)) {
        suanpan_error("create_new_mpc() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_mpc() needs a valid magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    vector<uword> node_tag, dof_tag;
    vector<double> weight_tag;
    while(!command.eof()) {
        double weight;
        uword dof, node;
        if(!get_input(command, node) || !get_input(command, dof) || !get_input(command, weight)) return SUANPAN_SUCCESS;
        node_tag.emplace_back(node);
        dof_tag.emplace_back(dof);
        weight_tag.emplace_back(weight);
    }

    domain->insert(make_shared<MPC>(tag, domain->get_current_step_tag(), amplitude, uvec(node_tag), uvec(dof_tag), vec(weight_tag), magnitude));

    return SUANPAN_SUCCESS;
}

int create_new_node(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned node_id;
    if(!get_input(command, node_id)) {
        suanpan_error("create_new_node() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    vector<double> coor;
    double X;
    while(get_input(command, X)) coor.push_back(X);

    if(!domain->insert(make_shared<Node>(node_id, vec(coor)))) suanpan_error("create_new_node() fails to insert Node %u.\n", node_id);

    return SUANPAN_SUCCESS;
}

int create_new_orientation(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string file_type;
    if(!get_input(command, file_type)) {
        suanpan_error("create_new_orientation() needs a valid type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_orientation() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    vec xyz(3);
    for(auto& I : xyz)
        if(!get_input(command, I)) {
            suanpan_error("create_new_orientation() needs a vector.\n");
            return SUANPAN_SUCCESS;
        }

    if(is_equal(file_type, "B3DL")) domain->insert(make_shared<B3DL>(tag, std::move(xyz)));
    else if(is_equal(file_type, "B3DC")) domain->insert(make_shared<B3DC>(tag, std::move(xyz)));

    return SUANPAN_SUCCESS;
}

int create_new_recorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_recorder() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    string file_type;
    if(!get_input(command, file_type)) {
        suanpan_error("create_new_recorder() needs a valid object type.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("create_new_recorder() needs a valid object type.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Eigen")) {
        if(!domain->insert(make_shared<EigenRecorder>(tag, is_equal(file_type[0], 'h')))) suanpan_error("create_new_recorder() fails to create a new eigen recorder.\n");
        return SUANPAN_SUCCESS;
    }

    string variable_type;
    if(!is_equal(object_type, "Amplitude") && !get_input(command, variable_type)) {
        suanpan_error("create_new_recorder() needs a valid recorder type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned interval = 1;

    while(true)
        if(const auto peek_value = command.peek(); is_equal(peek_value, '\t') || is_equal(peek_value, ' ')) command.ignore();
        else break;

    if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
        string tmp_string;
        get_input(command, tmp_string);
        if(!get_input(command, interval)) return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Frame")) {
        if(!domain->insert(make_shared<FrameRecorder>(tag, to_list(variable_type.c_str()), interval))) suanpan_error("create_new_recorder() fails to create a new frame recorder.\n");
        return SUANPAN_SUCCESS;
    }
    if(is_equal(object_type, "Visualisation")) {
        unsigned width = 6;
        if(!command.eof() && !get_input(command, width)) width = 6;
        if(!domain->insert(make_shared<VisualisationRecorder>(tag, to_list(variable_type.c_str()), interval, width))) suanpan_error("create_new_recorder() fails to create a new visualisation recorder.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned s_object_tag;
    vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(const auto use_hdf5 = is_equal(file_type[0], 'h'); is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5))) suanpan_error("create_new_recorder() fails to create a new node recorder.\n");
    else if(is_equal(object_type, "GroupNode") && !domain->insert(make_shared<GroupNodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5))) suanpan_error("create_new_recorder() fails to create a new group node recorder.\n");
    else if(is_equal(object_type, "Sum") && !domain->insert(make_shared<SumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5))) suanpan_error("create_new_recorder() fails to create a new summation recorder.\n");
    else if(is_equal(object_type, "GroupSum") && !domain->insert(make_shared<GroupSumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5))) suanpan_error("create_new_recorder() fails to create a new group summation recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5))) suanpan_error("create_new_recorder() fails to create a new element recorder.\n");
    else if(is_equal(object_type, "GroupElement") && !domain->insert(make_shared<GroupElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5))) suanpan_error("create_new_recorder() fails to create a new group element recorder.\n");
    else if(is_equal(object_type, "Amplitude") && !domain->insert(make_shared<AmplitudeRecorder>(tag, uvec(object_tag), OutputType::AMP, interval, true, use_hdf5))) suanpan_error("create_new_recorder() fails to create a new amplitude recorder.\n");
    else if(is_equal(object_type, "Global")) {
        bool flag;
        if(OutputType::K == to_list(variable_type.c_str())) flag = domain->insert(make_shared<GlobalStiffnessRecorder>(tag, interval, true, use_hdf5));
        else if(OutputType::M == to_list(variable_type.c_str())) flag = domain->insert(make_shared<GlobalMassRecorder>(tag, interval, true, use_hdf5));
        else flag = domain->insert(make_shared<GlobalRecorder>(tag, to_list(variable_type.c_str()), interval, true, use_hdf5));
        if(!flag) suanpan_error("create_new_hdf5recorder() fails to create a new global recorder.\n");
    }

    return SUANPAN_SUCCESS;
}

int create_new_plainrecorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_plainrecorder() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("create_new_plainrecorder() needs a valid object type.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Eigen")) {
        if(!domain->insert(make_shared<EigenRecorder>(tag, false))) suanpan_error("create_new_plainrecorder() fails to create a new eigen recorder.\n");
        return SUANPAN_SUCCESS;
    }

    string variable_type;
    if(!is_equal(object_type, "Amplitude") && !get_input(command, variable_type)) {
        suanpan_error("create_new_plainrecorder() needs a valid recorder type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned interval = 1;

    while(true)
        if(const auto peek_value = command.peek(); is_equal(peek_value, '\t') || is_equal(peek_value, ' ')) command.ignore();
        else break;

    if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
        string tmp_string;
        get_input(command, tmp_string);
        if(!get_input(command, interval)) return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Visualisation")) {
        unsigned width = 6;
        if(!command.eof() && !get_input(command, width)) width = 6;
        if(!domain->insert(make_shared<VisualisationRecorder>(tag, to_list(variable_type.c_str()), interval, width))) suanpan_error("create_new_recorder() fails to create a new visualisation recorder.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned s_object_tag;
    vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false))) suanpan_error("create_new_plainrecorder() fails to create a new node recorder.\n");
    else if(is_equal(object_type, "GroupNode") && !domain->insert(make_shared<GroupNodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false))) suanpan_error("create_new_plainrecorder() fails to create a new group node recorder.\n");
    else if(is_equal(object_type, "Sum") && !domain->insert(make_shared<SumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false))) suanpan_error("create_new_plainrecorder() fails to create a new summation recorder.\n");
    else if(is_equal(object_type, "GroupSum") && !domain->insert(make_shared<GroupSumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false))) suanpan_error("create_new_plainrecorder() fails to create a new group summation recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false))) suanpan_error("create_new_plainrecorder() fails to create a new element recorder.\n");
    else if(is_equal(object_type, "GroupElement") && !domain->insert(make_shared<GroupElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false))) suanpan_error("create_new_plainrecorder() fails to create a new group element recorder.\n");
    else if(is_equal(object_type, "Amplitude") && !domain->insert(make_shared<AmplitudeRecorder>(tag, uvec(object_tag), OutputType::AMP, interval, true, false))) suanpan_error("create_new_plainrecorder() fails to create a new amplitude recorder.\n");
    else if(is_equal(object_type, "Global")) {
        bool flag;
        if(OutputType::K == to_list(variable_type.c_str())) flag = domain->insert(make_shared<GlobalStiffnessRecorder>(tag, interval, true, false));
        else if(OutputType::M == to_list(variable_type.c_str())) flag = domain->insert(make_shared<GlobalMassRecorder>(tag, interval, true, false));
        else flag = domain->insert(make_shared<GlobalRecorder>(tag, to_list(variable_type.c_str()), interval, true, false));
        if(!flag) suanpan_error("create_new_hdf5recorder() fails to create a new global recorder.\n");
    }

    return SUANPAN_SUCCESS;
}

int create_new_hdf5recorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_hdf5recorder() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("create_new_hdf5recorder() needs a valid object type.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Eigen")) {
        if(!domain->insert(make_shared<EigenRecorder>(tag, true))) suanpan_error("create_new_hdf5recorder() fails to create a new eigen recorder.\n");
        return SUANPAN_SUCCESS;
    }

    string variable_type;
    if(!is_equal(object_type, "Amplitude") && !get_input(command, variable_type)) {
        suanpan_error("create_new_hdf5recorder() needs a valid recorder type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned interval = 1;

    while(true)
        if(const auto peek_value = command.peek(); is_equal(peek_value, '\t') || is_equal(peek_value, ' ')) command.ignore();
        else break;

    if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
        string tmp_string;
        get_input(command, tmp_string);
        if(!get_input(command, interval)) return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Frame")) {
        if(!domain->insert(make_shared<FrameRecorder>(tag, to_list(variable_type.c_str()), interval))) suanpan_error("create_new_recorder() fails to create a new frame recorder.\n");
        return SUANPAN_SUCCESS;
    }
    if(is_equal(object_type, "Visualisation")) {
        string para;
        unsigned width = 6;
        auto scale = 1.;
        while(!command.eof() && get_input(command, para))
            if(is_equal(para, "Width")) {
                if(!get_input(command, width)) {
                    width = 6;
                    suanpan_error("create_new_recorder() needs a proper width.\n");
                }
            }
            else if(is_equal(para, "Scale")) {
                if(!get_input(command, scale)) {
                    scale = 1.;
                    suanpan_error("create_new_recorder() needs a proper scale.\n");
                }
            }
        if(!domain->insert(make_shared<VisualisationRecorder>(tag, to_list(variable_type.c_str()), interval, width, scale))) suanpan_error("create_new_recorder() fails to create a new visualisation recorder.\n");
        return SUANPAN_SUCCESS;
    }

    uword s_object_tag = 0;
    vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true))) suanpan_error("create_new_hdf5recorder() fails to create a new node recorder.\n");
    else if(is_equal(object_type, "GroupNode") && !domain->insert(make_shared<GroupNodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true))) suanpan_error("create_new_hdf5recorder() fails to create a new group node recorder.\n");
    else if(is_equal(object_type, "Sum") && !domain->insert(make_shared<SumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true))) suanpan_error("create_new_hdf5recorder() fails to create a new summation recorder.\n");
    else if(is_equal(object_type, "GroupSum") && !domain->insert(make_shared<GroupSumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true))) suanpan_error("create_new_hdf5recorder() fails to create a new group summation recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true))) suanpan_error("create_new_hdf5recorder() fails to create a new element recorder.\n");
    else if(is_equal(object_type, "GroupElement") && !domain->insert(make_shared<GroupElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true))) suanpan_error("create_new_hdf5recorder() fails to create a new group element recorder.\n");
    else if(is_equal(object_type, "Amplitude") && !domain->insert(make_shared<AmplitudeRecorder>(tag, uvec(object_tag), OutputType::AMP, interval, true, true))) suanpan_error("create_new_hdf5recorder() fails to create a new amplitude recorder.\n");
    else if(is_equal(object_type, "Global")) {
        bool flag;
        if(OutputType::K == to_list(variable_type.c_str())) flag = domain->insert(make_shared<GlobalStiffnessRecorder>(tag, interval, true, true));
        else if(OutputType::M == to_list(variable_type.c_str())) flag = domain->insert(make_shared<GlobalMassRecorder>(tag, interval, true, true));
        else flag = domain->insert(make_shared<GlobalRecorder>(tag, to_list(variable_type.c_str()), interval, true, true));
        if(!flag) suanpan_error("create_new_hdf5recorder() fails to create a new global recorder.\n");
    }

    return SUANPAN_SUCCESS;
}

int create_new_rigidwall(const shared_ptr<DomainBase>& domain, istringstream& command, const bool finite, const bool penalty) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_rigidwall() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    vec origin(3), norm(3), edge(3);
    for(auto& I : origin) if(!get_input(command, I)) return SUANPAN_SUCCESS;
    for(auto& I : norm) if(!get_input(command, I)) return SUANPAN_SUCCESS;

    if(finite) for(auto& I : edge) if(!get_input(command, I)) return SUANPAN_SUCCESS;

    auto alpha = 1.;
    if(!command.eof() && !get_input(command, alpha)) {
        suanpan_error("create_new_rigidwall() needs a valid multiplier.\n");
        return SUANPAN_SUCCESS;
    }

    domain->insert(finite ? penalty ? make_shared<RigidWallPenalty>(tag, domain->get_current_step_tag(), 0, std::move(origin), std::move(norm), std::move(edge), alpha) : make_shared<RigidWallMultiplier>(tag, domain->get_current_step_tag(), 0, std::move(origin), std::move(norm), std::move(edge), alpha) : penalty ? make_shared<RigidWallPenalty>(tag, domain->get_current_step_tag(), 0, std::move(origin), std::move(norm), alpha) : make_shared<RigidWallMultiplier>(tag, domain->get_current_step_tag(), 0, std::move(origin), std::move(norm), alpha));

    return SUANPAN_SUCCESS;
}

int create_new_particlecollision2d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_particlecollision2d() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    auto space = 1.;
    if(!command.eof() && !get_input(command, space)) {
        suanpan_error("create_new_particlecollision2d() needs a valid spacing.\n");
        return SUANPAN_SUCCESS;
    }

    auto alpha = 1.;
    if(!command.eof() && !get_input(command, alpha)) {
        suanpan_error("create_new_particlecollision2d() needs a valid multiplier.\n");
        return SUANPAN_SUCCESS;
    }

    domain->insert(make_shared<ParticleCollision2D>(tag, domain->get_current_step_tag(), space, alpha));

    return SUANPAN_SUCCESS;
}

int create_new_particlecollision3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_particlecollision3d() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    auto space = 1.;
    if(!command.eof() && !get_input(command, space)) {
        suanpan_error("create_new_particlecollision3d() needs a valid spacing.\n");
        return SUANPAN_SUCCESS;
    }

    auto alpha = 1.;
    if(!command.eof() && !get_input(command, alpha)) {
        suanpan_error("create_new_particlecollision3d() needs a valid multiplier.\n");
        return SUANPAN_SUCCESS;
    }

    domain->insert(make_shared<ParticleCollision3D>(tag, domain->get_current_step_tag(), space, alpha));

    return SUANPAN_SUCCESS;
}

int create_new_solver(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string solver_type;
    if(!get_input(command, solver_type)) {
        suanpan_error("create_new_solver() requires solver type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_solver() requires a tag.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    if(is_equal(solver_type, "Newton")) { if(domain->insert(make_shared<Newton>(tag))) code = 1; }
    else if(is_equal(solver_type, "modifiedNewton") || is_equal(solver_type, "mNewton")) { if(domain->insert(make_shared<Newton>(tag, true))) code = 1; }
    else if(is_equal(solver_type, "BFGS")) { if(domain->insert(make_shared<BFGS>(tag))) code = 1; }
    else if(is_equal(solver_type, "LBFGS")) {
        auto max_history = 20;
        if(!command.eof() && !get_input(command, max_history)) {
            suanpan_error("create_new_solver() requires a valid maximum step for LBFGS algorithm.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<BFGS>(tag, max_history))) code = 1;
    }
    else if(is_equal(solver_type, "Ramm")) {
        auto arc_length = .1;
        string fixed_arc_length = "False";

        if(!command.eof() && !get_input(command, arc_length)) {
            suanpan_error("create_new_solver() requires a valid arc length.\n");
            return SUANPAN_SUCCESS;
        }
        if(!command.eof() && !get_input(command, fixed_arc_length)) {
            suanpan_error("create_new_solver() requires a valid arc length switch.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<Ramm>(tag, arc_length, is_true(fixed_arc_length)))) code = 1;
    }
    else if(is_equal(solver_type, "FEAST")) {
        unsigned eigen_number;
        if(!get_input(command, eigen_number)) {
            suanpan_error("create_new_solver() requires a valid number of frequencies.\n");
            return SUANPAN_SUCCESS;
        }

        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("create_new_solver() requires a valid radius.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<FEAST>(tag, eigen_number, radius, false))) code = 1;
    }
    else if(is_equal(solver_type, "QuadraticFEAST")) {
        unsigned eigen_number;
        if(!get_input(command, eigen_number)) {
            suanpan_error("create_new_solver() requires a valid number of frequencies.\n");
            return SUANPAN_SUCCESS;
        }

        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("create_new_solver() requires a valid radius.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<FEAST>(tag, eigen_number, radius, true))) code = 1;
    }
    else if(is_equal(solver_type, "DisplacementControl") || is_equal(solver_type, "MPDC")) { if(domain->insert(make_shared<MPDC>(tag))) code = 1; }
    else suanpan_error("create_new_solver() cannot identify solver type.\n");

    if(1 == code) {
        if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_solver_tag(tag);
        domain->set_current_solver_tag(tag);
    }
    else suanpan_error("create_new_solver() cannot create the new solver.\n");

    return SUANPAN_SUCCESS;
}

int create_new_step(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string step_type;
    if(!get_input(command, step_type)) {
        suanpan_error("create_new_step() requires step type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_step() requires a tag.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(step_type, "Frequency")) {
        auto eigen_number = 1;
        if(!command.eof() && !get_input(command, eigen_number)) {
            suanpan_error("create_new_step() reads a wrong number of eigenvalues.\n");
            return SUANPAN_SUCCESS;
        }
        if(domain->insert(make_shared<Frequency>(tag, eigen_number))) domain->set_current_step_tag(tag);
        else suanpan_error("create_new_step() cannot create the new step.\n");
    }
    else if(is_equal(step_type, "Buckling") || is_equal(step_type, "Buckle")) {
        if(domain->insert(make_shared<Buckle>(tag))) domain->set_current_step_tag(tag);
        else suanpan_error("create_new_step() cannot create the new step.\n");
    }
    else if(is_equal(step_type, "Optimization") || is_equal(step_type, "Optimisation")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            suanpan_error("create_new_step() reads a wrong time period.\n");
            return SUANPAN_SUCCESS;
        }
        if(domain->insert(make_shared<Optimization>(tag, time))) domain->set_current_step_tag(tag);
        else suanpan_error("create_new_step() cannot create the new step.\n");
    }
    else if(is_equal(step_type, "Static")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            suanpan_error("create_new_step() reads a wrong time period.\n");
            return SUANPAN_SUCCESS;
        }
        if(domain->insert(make_shared<Static>(tag, time))) domain->set_current_step_tag(tag);
        else suanpan_error("create_new_step() cannot create the new step.\n");
    }
    else if(is_equal(step_type, "Dynamic")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            suanpan_error("create_new_step() reads a wrong time period.\n");
            return SUANPAN_SUCCESS;
        }
        if(domain->insert(make_shared<Dynamic>(tag, time))) domain->set_current_step_tag(tag);
        else suanpan_error("create_new_step() cannot create the new step.\n");
    }
    else if(is_equal(step_type, "ArcLength")) {
        unsigned node;
        if(!get_input(command, node)) {
            suanpan_error("create_new_step() requires a node.\n");
            return SUANPAN_SUCCESS;
        }

        unsigned dof;
        if(!get_input(command, dof)) {
            suanpan_error("create_new_step() requires a dof.\n");
            return SUANPAN_SUCCESS;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("create_new_step() requires a magnitude.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<ArcLength>(tag, node, dof, magnitude))) domain->set_current_step_tag(tag);
        else suanpan_error("create_new_step() cannot create the new step.\n");
    }
    else suanpan_error("create_new_step() cannot identify step type.\n");

    return SUANPAN_SUCCESS;
}

int create_new_supportmotion(const shared_ptr<DomainBase>& domain, istringstream& command, const unsigned flag) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_supportmotion() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_supportmotion() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_supportmotion() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_supportmotion() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(const auto step_tag = domain->get_current_step_tag(); 0 == flag) { if(!domain->insert(make_shared<SupportDisplacement>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_supportmotion() fails to create new load.\n"); }
    else if(1 == flag) { if(!domain->insert(make_shared<SupportVelocity>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_supportmotion() fails to create new load.\n"); }
    else { if(!domain->insert(make_shared<SupportAcceleration>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_supportmotion() fails to create new load.\n"); }

    return SUANPAN_SUCCESS;
}

int test_material1d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material1d() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    double incre;
    if(!get_input(command, incre)) {
        suanpan_error("test_material1d() needs a valid step size.\n");
        return SUANPAN_SUCCESS;
    }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? static_cast<unsigned>(-step) : static_cast<unsigned>(step));

    if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester(material_proto->get_copy(), load_step, {incre});

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

    if(std::ofstream gnuplot("RESULT.plt"); gnuplot.is_open()) {
        gnuplot << "reset\n";
        gnuplot << "set term tikz size 14cm,10cm\n";
        gnuplot << "set output \"RESULT.tex\"\n";
        gnuplot << "unset key\n";
        gnuplot << "set xrange [*:*]\n";
        gnuplot << "set yrange [*:*]\n";
        gnuplot << "set xlabel \"input\"\n";
        gnuplot << "set ylabel \"output\"\n";
        gnuplot << "set grid\n";
        gnuplot << "plot \"RESULT.txt\" u 1:2 w l lw 2\n";
        gnuplot << "set output\n";

        gnuplot.close();
    }

    return SUANPAN_SUCCESS;
}

int test_material2d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material2d() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(3);
    for(auto& I : incre) {
        if(!get_input(command, I)) {
            suanpan_error("test_material2d() needs a valid step size.\n");
            return SUANPAN_SUCCESS;
        }
    }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? static_cast<unsigned>(-step) : static_cast<unsigned>(step));

    if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester(material_proto->get_copy(), load_step, incre);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

    return SUANPAN_SUCCESS;
}

int test_material3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material3d() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(6);
    for(auto& I : incre)
        if(!get_input(command, I)) {
            suanpan_error("test_material3d() needs a valid step size.\n");
            return SUANPAN_SUCCESS;
        }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? static_cast<unsigned>(-step) : static_cast<unsigned>(step));

    if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester(material_proto->get_copy(), load_step, incre);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_with_base3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material3dwithbase() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    vec base(6);
    for(auto& I : base)
        if(!get_input(command, I)) {
            suanpan_error("test_material3dwithbase() needs a valid step size.\n");
            return SUANPAN_SUCCESS;
        }

    vec incre(6);
    for(auto& I : incre)
        if(!get_input(command, I)) {
            suanpan_error("test_material3dwithbase() needs a valid step size.\n");
            return SUANPAN_SUCCESS;
        }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? static_cast<unsigned>(-step) : static_cast<unsigned>(step));

    if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester(material_proto->get_copy(), load_step, incre, base);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_by_load1d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material1d() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    double incre;
    if(!get_input(command, incre)) {
        suanpan_error("test_material1d() needs a valid step size.\n");
        return SUANPAN_SUCCESS;
    }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? static_cast<unsigned>(-step) : static_cast<unsigned>(step));

    if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester_by_load(material_proto->get_copy(), load_step, {incre});

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

    if(std::ofstream gnuplot("RESULT.plt"); gnuplot.is_open()) {
        gnuplot << "reset\n";
        gnuplot << "set term tikz size 14cm,10cm\n";
        gnuplot << "set output \"RESULT.tex\"\n";
        gnuplot << "unset key\n";
        gnuplot << "set xrange [*:*]\n";
        gnuplot << "set yrange [*:*]\n";
        gnuplot << "set xlabel \"input\"\n";
        gnuplot << "set ylabel \"output\"\n";
        gnuplot << "set grid\n";
        gnuplot << "plot \"RESULT.txt\" u 1:2 w l lw 2\n";
        gnuplot << "set output\n";
    }

    return SUANPAN_SUCCESS;
}

int test_material_by_load2d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material2d() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(3);
    for(auto& I : incre)
        if(!get_input(command, I)) {
            suanpan_error("test_material2d() needs a valid step size.\n");
            return SUANPAN_SUCCESS;
        }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? static_cast<unsigned>(-step) : static_cast<unsigned>(step));

    if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester_by_load(material_proto->get_copy(), load_step, incre);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_by_load3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material3d() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(6);
    for(auto& I : incre)
        if(!get_input(command, I)) {
            suanpan_error("test_material3d() needs a valid step size.\n");
            return SUANPAN_SUCCESS;
        }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? static_cast<unsigned>(-step) : static_cast<unsigned>(step));

    if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester_by_load(material_proto->get_copy(), load_step, incre);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_by_load_with_base3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material3dwithbase() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    vec base(6);
    for(auto& I : base)
        if(!get_input(command, I)) {
            suanpan_error("test_material3dwithbase() needs a valid step size.\n");
            return SUANPAN_SUCCESS;
        }

    vec incre(6);
    for(auto& I : incre)
        if(!get_input(command, I)) {
            suanpan_error("test_material3dwithbase() needs a valid step size.\n");
            return SUANPAN_SUCCESS;
        }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? static_cast<unsigned>(-step) : static_cast<unsigned>(step));

    if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester_by_load(material_proto->get_copy(), load_step, incre, base);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_by_strain_history(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material_by_strain_history() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    string history_file;
    if(!get_input(command, history_file)) {
        suanpan_error("test_material_by_strain_history() needs a valid history file name.\n");
        return SUANPAN_SUCCESS;
    }

    mat strain_history;
    if(!strain_history.load(history_file) || !domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester_by_strain_history(material_proto->get_copy(), strain_history);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#else
    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");
#endif

    return SUANPAN_SUCCESS;
}

int test_material_by_stress_history(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material_by_stress_history() needs a valid material tag.\n");
        return SUANPAN_SUCCESS;
    }

    string history_file;
    if(!get_input(command, history_file)) {
        suanpan_error("test_material_by_stress_history() needs a valid history file name.\n");
        return SUANPAN_SUCCESS;
    }

    mat stress_history;
    if(!stress_history.load(history_file) || !domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester_by_stress_history(material_proto->get_copy(), stress_history);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#else
    if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");
#endif

    return SUANPAN_SUCCESS;
}

int set_property(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string property_id;
    if(!get_input(command, property_id)) {
        suanpan_error("set_property() need a property type.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(property_id, "output_folder")) {
        string value;

        if(!get_input(command, value)) {
            suanpan_error("set_property() need a valid value.\n");
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
                    suanpan_error("cannot create %s.\n", new_path.generic_string().c_str());
                    return SUANPAN_SUCCESS;
                }
            }

            SUANPAN_OUTPUT = canonical(new_path);
        }

        suanpan_info("%s\n", SUANPAN_OUTPUT.generic_string().c_str());
        return SUANPAN_SUCCESS;
    }
    if(is_equal(property_id, "num_threads")) {
        if(int value; get_input(command, value)) SUANPAN_NUM_THREADS = value;
        else suanpan_error("set_property() need a valid value.\n");

        return SUANPAN_SUCCESS;
    }
    if(is_equal(property_id, "screen_output")) {
        if(string value; get_input(command, value)) SUANPAN_PRINT = is_true(value);
        else suanpan_error("set_property() need a valid value.\n");

        return SUANPAN_SUCCESS;
    }

    if(domain->get_current_step_tag() == 0) return SUANPAN_SUCCESS;

    const auto& t_step = domain->get_current_step();

    if(is_equal(property_id, "color_model")) {
        if(string value; !get_input(command, value)) suanpan_error("set_property() need a valid value.\n");
        else if(is_equal("WP", value)) domain->set_color_model(ColorMethod::WP);
        else if(is_equal("MIS", value)) domain->set_color_model(ColorMethod::MIS);
        else domain->set_color_model(ColorMethod::OFF);
    }
    else if(is_equal(property_id, "fixed_step_size")) {
        string value;
        get_input(command, value) ? t_step->set_fixed_step_size(is_true(value)) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "symm_mat")) {
        string value;
        get_input(command, value) ? t_step->set_symm(is_true(value)) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "band_mat")) {
        string value;
        get_input(command, value) ? t_step->set_band(is_true(value)) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "sparse_mat")) {
        string value;
        get_input(command, value) ? t_step->set_sparse(is_true(value)) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "iterative_refinement")) {
        if(unsigned value; get_input(command, value)) t_step->set_refinement(value);
        else suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "system_solver")) {
        if(string value; !get_input(command, value)) suanpan_error("set_property() need a valid value.\n");
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
        else suanpan_error("set_property() need a valid solver id.\n");
    }
    else if(is_equal(property_id, "precision")) {
        if(string value; !get_input(command, value)) suanpan_error("set_property() need a valid value.\n");
        else if(is_equal(value, "DOUBLE") || is_equal(value, "FULL")) t_step->set_precision(Precision::FULL);
        else if(is_equal(value, "SINGLE") || is_equal(value, "MIXED")) t_step->set_precision(Precision::MIXED);
        else suanpan_error("set_property() need a valid precision.\n");
    }
    else if(is_equal(property_id, "tolerance")) {
        double value;
        get_input(command, value) ? t_step->set_tolerance(value) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "ini_step_size")) {
        double step_time;
        get_input(command, step_time) ? t_step->set_ini_step_size(step_time) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "min_step_size")) {
        double step_time;
        get_input(command, step_time) ? t_step->set_min_step_size(step_time) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "max_step_size")) {
        double step_time;
        get_input(command, step_time) ? t_step->set_max_step_size(step_time) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "max_iteration")) {
        unsigned max_number;
        get_input(command, max_number) ? t_step->set_max_substep(max_number) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "eigen_number")) {
        if(unsigned eigen_number; get_input(command, eigen_number)) {
            if(const auto eigen_step = std::dynamic_pointer_cast<Frequency>(t_step); nullptr == eigen_step) suanpan_error("set_property() cannot set eigen number for noneigen step.\n");
            else eigen_step->set_eigen_number(eigen_number);
        }
        else suanpan_error("set_property() need a valid eigen number.\n");
    }
    else if(is_equal(property_id, "constraint_multiplier")) {
        double value;
        get_input(command, value) ? set_constraint_multiplier(value) : suanpan_error("set_property() need a valid value.\n");
    }
    else if(is_equal(property_id, "load_multiplier")) {
        double value;
        get_input(command, value) ? set_load_multiplier(value) : suanpan_error("set_property() need a valid value.\n");
    }
#ifdef SUANPAN_MKL
    else if(is_equal(property_id, "fgmres_tolerance")) {
        double value;
        get_input(command, value) ? set_fgmres_tolerance(value) : suanpan_error("set_property() need a valid value.\n");
    }
#endif

    return SUANPAN_SUCCESS;
}

int print_info(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string object_type;
    if(!get_input(command, object_type)) {
        suanpan_error("print_info() needs object type.\n");
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
    else if(is_equal(object_type, "eigenvalue")) {
        domain->get_factory()->get_eigenvalue().print("Eigenvalues:");
        suanpan_info("\n");
    }
    else if(is_equal(object_type, "output_folder")) suanpan_info("%s\n", SUANPAN_OUTPUT.generic_string().c_str());

    return SUANPAN_SUCCESS;
}

int run_example() {
    const auto new_model = make_shared<Bead>();

    suanpan_info("====================================================\n");
    suanpan_info("-> A Minimum Example: Elastic Truss Under Tension <-\n");
    suanpan_info("====================================================\n");

    constexpr auto wait_time = 1000;

    auto run_command = [&](const string& command_string) {
        suanpan_info("\t%s\n", command_string.c_str());
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

int print_command() {
    suanpan_info("The available commands are listed. Please check online manual for reference. https://tlcfem.gitbook.io/suanpan-manual/\n");

    constexpr auto format = "    %-30s  %s\n";
    suanpan_info(format, "acceleration", "define acceleration");
    suanpan_info(format, "amplitude", "define amplitude");
    suanpan_info(format, "analyze/analyse", "analyse the model");
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
    suanpan_info(format, "exit/quit", "exit program");
    suanpan_info(format, "file/load", "load external files");
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
    suanpan_info(format, "version", "print version information");

    return SUANPAN_SUCCESS;
}

int execute_command(istringstream& command) {
#ifdef SUANPAN_MSVC
    std::wstringstream terminal_command;
    terminal_command << command.str().substr(command.tellg()).c_str();
    return _wsystem(terminal_command.str().c_str());
#else
    std::stringstream terminal_command;
    terminal_command << command.str().substr(command.tellg()).c_str();
    return system(terminal_command.str().c_str());
#endif
}
