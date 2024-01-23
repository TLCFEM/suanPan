/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "SingleSection2D.h"
#include <Domain/DomainBase.h>
#include <Section/Section.h>

SingleSection2D::SingleSection2D(const unsigned T, const unsigned NT, const unsigned ST)
    : SectionElement2D(T, s_node, s_dof, uvec{NT}, uvec{ST}, false, {}) {}

int SingleSection2D::initialize(const shared_ptr<DomainBase>& D) {
    s_section = suanpan::make_copy(D->get<Section>(section_tag(0)));

    initial_stiffness = s_section->get_initial_stiffness();

    return SUANPAN_SUCCESS;
}

int SingleSection2D::update_status() {
    if(SUANPAN_SUCCESS != s_section->update_trial_status(get_trial_displacement())) return SUANPAN_FAIL;

    trial_stiffness = s_section->get_trial_stiffness();

    trial_resistance = s_section->get_trial_resistance();

    return SUANPAN_SUCCESS;
}

int SingleSection2D::commit_status() { return s_section->commit_status(); }

int SingleSection2D::clear_status() { return s_section->clear_status(); }

int SingleSection2D::reset_status() { return s_section->reset_status(); }

void SingleSection2D::print() {
    suanpan_info("A SingleSection2D element that represents a section which can be used for section analysis.\n");
    if(s_section) s_section->print();
}
