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

#include "ElementExample.h"
#include <Domain/DOF.h>
#include <Domain/DomainBase.h>
#include <Material/Material.h>
#include <Toolbox/utility.h>

SUANPAN_EXPORT void new_elementexample(unique_ptr<Element>& return_obj, std::istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    std::vector<uword> node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("Three valid nodes are required.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_debug("Unit thickness assumed.\n");
    else if(!get_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    return_obj = make_unique<ElementExample>(tag, node_tag, material_tag, thickness);
}

ElementExample::ElementExample(const unsigned T, uvec&& NT, const unsigned MT, const double TH)
    : Element(T, m_node, m_dof, std::move(NT), uvec{MT}, false, MaterialType::D2, {DOF::U1, DOF::U2})
    , thickness(TH) {}

int ElementExample::initialize(const shared_ptr<DomainBase>& D) {
    m_material = D->get<Material>(material_tag(0))->get_copy();

    mat ele_coor(m_node, m_node, fill::ones);
    ele_coor.cols(1, 2) = get_coordinate(2);

    access::rw(area) = .5 * det(ele_coor);

    const mat inv_coor = inv(ele_coor);

    strain_mat.zeros(3, m_size);
    for(unsigned i = 0, j = 0, k = 1; i < 3; ++i, j += m_dof, k += m_dof) {
        strain_mat(2, k) = strain_mat(0, j) = inv_coor(1, i);
        strain_mat(2, j) = strain_mat(1, k) = inv_coor(2, i);
    }

    trial_stiffness = current_stiffness = initial_stiffness = strain_mat.t() * m_material->get_initial_stiffness() * strain_mat * area * thickness;

    initial_mass.zeros(m_size, m_size);
    const rowvec n = mean(ele_coor) * inv_coor;
    const mat t_mass = n.t() * n * area * thickness * m_material->get_density();
    initial_mass(uvec{1, 3, 5}, uvec{1, 3, 5}) = t_mass;
    initial_mass(uvec{0, 2, 4}, uvec{0, 2, 4}) = t_mass;

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

int ElementExample::update_status() {
    if(m_material->update_trial_status(strain_mat * get_trial_displacement()) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stiffness = strain_mat.t() * m_material->get_trial_stiffness() * strain_mat * area * thickness;
    trial_resistance = strain_mat.t() * m_material->get_trial_stress() * area * thickness;

    return SUANPAN_SUCCESS;
}

int ElementExample::commit_status() { return m_material->commit_status(); }

int ElementExample::clear_status() { return m_material->clear_status(); }

int ElementExample::reset_status() { return m_material->reset_status(); }

void ElementExample::print() {
    suanpan_info("An element example based on CP3 element using the following material model.\n");
    if(m_material) m_material->print();
}
