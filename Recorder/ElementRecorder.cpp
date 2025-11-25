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

#include "ElementRecorder.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Element/Element.h>

void ElementRecorder::initialize(const shared_ptr<DomainBase>& D) {
    update_tag(D);

    std::vector<uword> pool;
    pool.reserve(object_tag.n_elem);
    for(const auto I : object_tag)
        if(!D->find<Element>(I) || !D->get<Element>(I)->is_active())
            suanpan_warning("Element {} is not available/active, removed from recorder {}.\n", I, get_tag());
        else pool.emplace_back(I);

    object_tag = pool;
    data_pool.resize(object_tag.n_elem);
}

void ElementRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    if(OutputType::K == variable_type) {
        for(auto I = 0llu; I < object_tag.n_elem; ++I)
            if(const auto& t_element = D->get<Element>(object_tag(I)); t_element->is_active() && t_element->is_local) insert({vectorise(t_element->get_current_stiffness())}, I);
    }
    else if(OutputType::M == variable_type) {
        for(auto I = 0llu; I < object_tag.n_elem; ++I)
            if(const auto& t_element = D->get<Element>(object_tag(I)); t_element->is_active() && t_element->is_local) insert({vectorise(t_element->get_current_mass())}, I);
    }
    else
        for(auto I = 0llu; I < object_tag.n_elem; ++I)
            if(const auto& t_element = D->get<Element>(object_tag(I)); t_element->is_active() && t_element->is_local) insert(t_element->record(variable_type), I);

    insert(D->get_factory()->get_current_time());
}

void ElementRecorder::print() { suanpan_info("An element recorder.\n"); }

const uvec& GroupElementRecorder::update_tag(const shared_ptr<DomainBase>& D) { return object_tag = D->flatten_group(reference_tag); }

void GroupElementRecorder::print() { suanpan_info("An element recorder based on groups.\n"); }
