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
    std::vector<uword> pool;
    pool.reserve(get_object_tag().n_elem);
    for(const auto I : get_object_tag())
        if(!D->find<Element>(I) || !D->get<Element>(I)->is_active())
            suanpan_warning("Element {} is not available/active, removed from recorder {}.\n", I, get_tag());
        else pool.emplace_back(I);

    set_object_tag(pool);

    access::rw(get_data_pool()).resize(get_object_tag().n_elem);
}

void ElementRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    auto& obj_tag = get_object_tag();

    if(OutputType::K == get_variable_type()) {
        for(unsigned I = 0; I < obj_tag.n_elem; ++I)
            if(const auto& t_element = D->get<Element>(obj_tag(I)); t_element->is_active() && t_element->is_local) insert({vectorise(t_element->get_current_stiffness())}, I);
    }
    else if(OutputType::M == get_variable_type()) {
        for(unsigned I = 0; I < obj_tag.n_elem; ++I)
            if(const auto& t_element = D->get<Element>(obj_tag(I)); t_element->is_active() && t_element->is_local) insert({vectorise(t_element->get_current_mass())}, I);
    }
    else
        for(unsigned I = 0; I < obj_tag.n_elem; ++I)
            if(const auto& t_element = D->get<Element>(obj_tag(I)); t_element->is_active() && t_element->is_local) insert(t_element->record(get_variable_type()), I);

    if(if_record_time()) insert(D->get_factory()->get_current_time());
}

void ElementRecorder::print() {
    suanpan_info("An element recorder.\n");
}
