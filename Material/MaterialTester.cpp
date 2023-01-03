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

#include "MaterialTester.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>
#include <Toolbox/utility.h>

bool initialise_material(const shared_ptr<Material>& obj, const uword size) {
    if(!obj->is_initialized()) {
        obj->initialize_base(nullptr);
        obj->initialize(nullptr);
        obj->set_initialized(true);
    }

    if(obj->get_material_type() == MaterialType::D1 && 1 != size) {
        suanpan_error("The tester cannot be applied to the given material model.\n");
        return false;
    }
    if(obj->get_material_type() == MaterialType::D2 && 3 != size) {
        suanpan_error("The tester cannot be applied to the given material model.\n");
        return false;
    }
    if(obj->get_material_type() == MaterialType::D3 && 6 != size) {
        suanpan_error("The tester cannot be applied to the given material model.\n");
        return false;
    }

    return true;
}

mat material_tester(const shared_ptr<Material>& obj, const std::vector<unsigned>& idx, const vec& incre) {
    if(!initialise_material(obj, incre.n_elem)) return {};

    unsigned total_size = 1;
    for(const auto& I : idx) total_size += I;

    mat response(total_size, 2 * incre.n_elem, fill::zeros);

    const span span_a(0, incre.n_elem - 1);
    const span span_b(incre.n_elem, 2 * incre.n_elem - 1);

    auto current_pos = 0;

    response(current_pos, span_a) = obj->get_current_strain().t();
    response(current_pos++, span_b) = obj->get_current_stress().t();

    auto flag = SUANPAN_SUCCESS;
    auto incre_strain = incre;
    for(const auto& I : idx) {
        for(unsigned J = 0; J < I; ++J) {
            if((flag = obj->update_incre_status(incre_strain)) != SUANPAN_SUCCESS) break;
            obj->commit_status();
            response(current_pos, span_a) = obj->get_current_strain().t();
            response(current_pos++, span_b) = obj->get_current_stress().t();
        }
        if(SUANPAN_SUCCESS != flag) break;
        incre_strain = -incre_strain;
    }

    obj->print();
    obj->reset_status();
    obj->clear_status();

    return response;
}

mat material_tester(const shared_ptr<Material>& obj, const std::vector<unsigned>& idx, const vec& incre, const vec& base) {
    if(!initialise_material(obj, incre.n_elem)) return {};

    unsigned total_size = 2;
    for(const auto& I : idx) total_size += I;

    mat response(total_size, 2 * incre.n_elem, fill::zeros);

    const span span_a(0, incre.n_elem - 1);
    const span span_b(incre.n_elem, 2 * incre.n_elem - 1);

    auto current_pos = 0;

    response(current_pos, span_a) = obj->get_current_strain().t();
    response(current_pos++, span_b) = obj->get_current_stress().t();

    if(obj->update_incre_status(base) != SUANPAN_SUCCESS) return response;

    obj->commit_status();

    response(current_pos, span_a) = obj->get_current_strain().t();
    response(current_pos++, span_b) = obj->get_current_stress().t();

    auto flag = SUANPAN_SUCCESS;
    auto incre_strain = incre;
    for(const auto& I : idx) {
        for(unsigned J = 0; J < I; ++J) {
            if((flag = obj->update_incre_status(incre_strain)) != SUANPAN_SUCCESS) break;
            obj->commit_status();
            response(current_pos, span_a) = obj->get_current_strain().t();
            response(current_pos++, span_b) = obj->get_current_stress().t();
        }
        if(SUANPAN_SUCCESS != flag) break;
        incre_strain = -incre_strain;
    }

    obj->print();
    obj->reset_status();
    obj->clear_status();

    return response;
}

mat material_tester_by_load(const shared_ptr<Material>& obj, const std::vector<unsigned>& idx, const vec& incre) {
    if(!initialise_material(obj, incre.n_elem)) return {};

    unsigned total_size = 1;
    for(const auto& I : idx) total_size += I;

    mat response(total_size, 2 * incre.n_elem, fill::zeros);

    const span span_a(0, incre.n_elem - 1);
    const span span_b(incre.n_elem, 2 * incre.n_elem - 1);

    auto current_pos = 0;

    response(current_pos, span_a) = obj->get_current_strain().t();
    response(current_pos++, span_b) = obj->get_current_stress().t();

    auto info = SUANPAN_SUCCESS;
    auto incre_load = incre;
    vec total_load = zeros(size(incre));
    for(const auto& I : idx) {
        for(unsigned J = 0; J < I; ++J) {
            total_load += incre_load;
            auto counter = 0;
            while(true) {
                const vec incre_strain = solve(obj->get_trial_stiffness(), total_load - obj->get_trial_stress());
                const auto error = norm(incre_strain);
                suanpan_debug("Local iteration error: {:.5E}.\n", error);
                if(error < 1E-12) break;
                if(++counter == 10 || obj->update_trial_status(obj->get_trial_strain() + incre_strain) != SUANPAN_SUCCESS) {
                    info = SUANPAN_FAIL;
                    break;
                }
            }
            if(SUANPAN_FAIL == info) break;
            obj->commit_status();
            response(current_pos, span_a) = obj->get_current_strain().t();
            response(current_pos++, span_b) = obj->get_current_stress().t();
        }
        if(SUANPAN_FAIL == info) break;
        incre_load = -incre_load;
    }

    obj->print();
    obj->reset_status();
    obj->clear_status();

    return response;
}

mat material_tester_by_load(const shared_ptr<Material>& obj, const std::vector<unsigned>& idx, const vec& incre, const vec& base) {
    if(!initialise_material(obj, incre.n_elem)) return {};

    unsigned total_size = 2;
    for(const auto& I : idx) total_size += I;

    mat response(total_size, 2 * incre.n_elem, fill::zeros);

    const span span_a(0, incre.n_elem - 1);
    const span span_b(incre.n_elem, 2 * incre.n_elem - 1);

    auto current_pos = 0;

    response(current_pos, span_a) = obj->get_current_strain().t();
    response(current_pos++, span_b) = obj->get_current_stress().t();

    auto incre_load = incre;
    auto total_load = base;
    auto info = SUANPAN_SUCCESS;
    auto counter = 0;
    while(true) {
        const vec incre_strain = solve(obj->get_trial_stiffness(), total_load - obj->get_trial_stress());
        const auto error = norm(incre_strain);
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < 1E-12) break;
        if(++counter == 10 || obj->update_trial_status(obj->get_trial_strain() + incre_strain) != SUANPAN_SUCCESS) {
            info = SUANPAN_FAIL;
            break;
        }
    }
    obj->commit_status();
    response(current_pos, span_a) = obj->get_current_strain().t();
    response(current_pos++, span_b) = obj->get_current_stress().t();

    if(SUANPAN_SUCCESS == info)
        for(const auto& I : idx) {
            for(unsigned J = 0; J < I; ++J) {
                total_load += incre_load;
                counter = 0;
                while(true) {
                    const vec incre_strain = solve(obj->get_trial_stiffness(), total_load - obj->get_trial_stress());
                    const auto error = norm(incre_strain);
                    suanpan_debug("Local iteration error: {:.5E}.\n", error);
                    if(error <= 1E-12) break;
                    if(++counter == 10 || obj->update_trial_status(obj->get_trial_strain() + incre_strain) != SUANPAN_SUCCESS) {
                        info = SUANPAN_FAIL;
                        break;
                    }
                }
                if(SUANPAN_FAIL == info) break;
                obj->commit_status();
                response(current_pos, span_a) = obj->get_current_strain().t();
                response(current_pos++, span_b) = obj->get_current_stress().t();
            }
            if(SUANPAN_FAIL == info) break;
            incre_load = -incre_load;
        }

    obj->print();
    obj->reset_status();
    obj->clear_status();

    return response;
}

mat material_tester_by_strain_history(const shared_ptr<Material>& obj, const mat& history) {
    if(!initialise_material(obj, history.n_cols)) return {};

    mat response(size(history));

    for(auto I = 0llu; I < history.n_rows; ++I) {
        if(SUANPAN_SUCCESS != obj->update_trial_status(history.row(I).t())) break;
        obj->commit_status();
        response.row(I) = obj->get_current_stress();
    }

    obj->print();
    obj->reset_status();
    obj->clear_status();

    return response;
}

mat material_tester_by_stress_history(const shared_ptr<Material>& obj, const mat& history) {
    if(!initialise_material(obj, history.n_cols)) return {};

    mat response(size(history));

    for(auto I = 0llu; I < history.n_rows; ++I) {
        auto counter = 0;
        auto flag = false;
        auto strain = obj->get_current_strain();
        while(true) {
            if(20 == ++counter) {
                flag = true;
                break;
            }
            const vec incre_strain = solve(obj->get_trial_stiffness(), history.row(I).t() - obj->get_trial_stress());
            const auto error = norm(incre_strain);
            suanpan_debug("Local iteration error: {:.5E}.\n", error);
            if(error <= 1E-12) break;
            strain += incre_strain;
            if(SUANPAN_SUCCESS != obj->update_trial_status(strain)) {
                flag = true;
                break;
            }
        }

        if(flag) break;

        obj->commit_status();
        response.row(I) = obj->get_current_strain();
    }

    obj->print();
    obj->reset_status();
    obj->clear_status();

    return response;
}

int test_material1d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    double incre;
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
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
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");

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
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(3);
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
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
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");

    return SUANPAN_SUCCESS;
}

int test_material3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(6);
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
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
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_with_base3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec base(6);
    if(!get_input(command, base)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(6);
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
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
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_by_load1d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    double incre;
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
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
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");

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
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(3);
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
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
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_by_load3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(6);
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
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
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_by_load_with_base3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec base(6);
    if(!get_input(command, base)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(6);
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
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
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#endif

    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");

    return SUANPAN_SUCCESS;
}

int test_material_by_strain_history(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    string history_file;
    if(!get_input(command, history_file)) {
        suanpan_error("A valid history file name is required.\n");
        return SUANPAN_SUCCESS;
    }

    mat strain_history;
    if(!strain_history.load(history_file, raw_ascii) || !domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester_by_strain_history(material_proto->get_copy(), strain_history);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#else
    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");
#endif

    return SUANPAN_SUCCESS;
}

int test_material_by_stress_history(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    string history_file;
    if(!get_input(command, history_file)) {
        suanpan_error("A valid history file name is required.\n");
        return SUANPAN_SUCCESS;
    }

    mat stress_history;
    if(!stress_history.load(history_file, raw_ascii) || !domain->find_material(material_tag)) return SUANPAN_SUCCESS;

    auto& material_proto = domain->get_material(material_tag);

    if(!material_proto->is_initialized()) {
        material_proto->initialize_base(domain);
        material_proto->initialize(domain);
        material_proto->set_initialized(true);
    }

    const auto result = material_tester_by_stress_history(material_proto->get_copy(), stress_history);

#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#else
    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");
#endif

    return SUANPAN_SUCCESS;
}
