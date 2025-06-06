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

#include "SectionTester.h"

#include <Domain/DomainBase.h>
#include <Section/Section.h>
#include <Toolbox/misc.h>
#include <Toolbox/utility.h>

bool initialise_section(const shared_ptr<DomainBase>& domain, const unique_ptr<Section>& obj, const uword size) {
    domain->initialize_material();
    if(!obj->is_initialized()) {
        if(SUANPAN_SUCCESS != obj->initialize_base(domain)) return false;
        if(SUANPAN_SUCCESS != obj->initialize(domain)) return false;
        obj->set_initialized(true);
    }

    if(static_cast<uword>(obj->get_section_type()) == size) return true;

    suanpan_error("The tester cannot be applied to the given section model.\n");
    return false;
}

mat section_tester(const unique_ptr<Section>& obj, const std::vector<unsigned>& idx, const vec& incre) {
    unsigned total_size = 1;
    for(const auto& I : idx) total_size += I;

    mat response(total_size, 2 * incre.n_elem, fill::zeros);

    const span span_a(0, incre.n_elem - 1);
    const span span_b(incre.n_elem, 2 * incre.n_elem - 1);

    auto current_pos = 0;

    response(current_pos, span_a) = obj->get_current_deformation().t();
    response(current_pos++, span_b) = obj->get_current_resistance().t();

    auto flag = SUANPAN_SUCCESS;
    auto incre_strain = incre;
    for(const auto& I : idx) {
        for(unsigned J = 0; J < I; ++J) {
            if((flag = obj->update_incre_status(incre_strain)) != SUANPAN_SUCCESS) break;
            obj->commit_status();
            response(current_pos, span_a) = obj->get_current_deformation().t();
            response(current_pos++, span_b) = obj->get_current_resistance().t();
        }
        if(SUANPAN_SUCCESS != flag) break;
        incre_strain = -incre_strain;
    }

    obj->print();
    obj->reset_status();
    obj->clear_status();

    return response;
}

mat section_tester_by_deformation_history(const unique_ptr<Section>& obj, const mat& history) {
    mat response(size(history));

    for(auto I = 0llu; I < history.n_rows; ++I) {
        if(SUANPAN_SUCCESS != obj->update_trial_status(history.row(I).t())) break;
        obj->commit_status();
        response.row(I) = obj->get_current_resistance().t();
    }

    obj->print();
    obj->reset_status();
    obj->clear_status();

    return response;
}

int test_section(const shared_ptr<DomainBase>& domain, std::istringstream& command, const unsigned size) {
    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("A valid section tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec incre(size);
    if(!get_input(command, incre)) {
        suanpan_error("A valid step size is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(static_cast<unsigned>(std::abs(step)));

    if(!domain->find_section(section_tag)) return SUANPAN_SUCCESS;

    const auto section = domain->get_section(section_tag)->get_copy();

    if(!initialise_section(domain, section, size)) return SUANPAN_SUCCESS;

    save_result(section_tester(section, load_step, incre));

    return SUANPAN_SUCCESS;
}

int test_section_by_deformation_history(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("A valid section tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::string history_file;
    if(!get_input(command, history_file)) {
        suanpan_error("A valid history file name is required.\n");
        return SUANPAN_SUCCESS;
    }

    mat deformation_history;
    if(!deformation_history.load(history_file, raw_ascii) || !domain->find_section(section_tag)) return SUANPAN_SUCCESS;

    const auto section = domain->get_section(section_tag)->get_copy();

    if(!initialise_section(domain, section, deformation_history.n_cols)) return SUANPAN_SUCCESS;

    save_result(section_tester_by_deformation_history(section, deformation_history));

    return SUANPAN_SUCCESS;
}
