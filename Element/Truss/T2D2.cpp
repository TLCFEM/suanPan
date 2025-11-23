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

#include "T2D2.h"

#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

T2D2::T2D2(const unsigned T, uvec&& N, const unsigned M, const double A, const bool F, const bool UA, const bool LS, const double FR)
    : MaterialElement1D(T, t_node, t_dof, std::move(N), uvec{M}, F, {Node::DOF::U1, Node::DOF::U2})
    , area(A)
    , rigidity(FR)
    , t_trans(F ? std::make_unique<T2DC>() : std::make_unique<T2DL>())
    , update_area(UA)
    , log_strain(LS) {}

int T2D2::initialize(const shared_ptr<DomainBase>& D) {
    t_trans->set_element_ptr(this);

    access::rw(length) = t_trans->get_length();

    if(rigidity > 0.) access::rw(euler) = -pow(datum::pi / length, 2.) * rigidity;

    t_material = D->get<Material>(material_tag(0))->get_copy();

    trial_stiffness = current_stiffness = initial_stiffness = t_trans->to_global_stiffness_mat(area / length * t_material->get_initial_stiffness());

    if(const auto t_density = t_material->get_density(); t_density > 0.) trial_mass = current_mass = initial_mass = t_trans->to_global_mass_mat(t_density * area);

    return SUANPAN_SUCCESS;
}

int T2D2::update_status() {
    auto new_area = area;

    t_trans->update_status();

    if(nlgeom) {
        const auto new_length = t_trans->get_length();

        double d_strain;

        if(t_material->update_trial_status(log_strain ? log((d_strain = new_length) / length) : new_length / (d_strain = length) - 1.) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        const auto t_stiff = update_area ? -t_material->get_trial_stress().at(0) * (new_area *= length / new_length) / new_length : 0.;

        trial_stiffness = t_trans->to_global_stiffness_mat(t_material->get_trial_stiffness() * new_area / d_strain + t_stiff);
        trial_geometry = t_trans->to_global_geometry_mat(new_area / new_length * t_material->get_trial_stress());
    }
    else {
        if(t_material->update_trial_status(t_trans->to_local_vec(get_trial_displacement()) / length) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness = t_trans->to_global_stiffness_mat(area / length * t_material->get_trial_stiffness());
    }

    const vec axial_force = new_area * t_material->get_trial_stress();
    trial_resistance = t_trans->to_global_vec(axial_force);

    if(euler < 0. && axial_force(0) <= euler) {
        suanpan_error("Element {} exceeds Euler buckling limit.\n", get_tag());
        return SUANPAN_FAIL;
    }

    suanpan_assert([&] { if(!trial_stiffness.is_finite() || !trial_resistance.is_finite()) throw std::invalid_argument("infinite number detected"); });

    return SUANPAN_SUCCESS;
}

int T2D2::commit_status() {
    t_trans->commit_status();

    return t_material->commit_status();
}

int T2D2::clear_status() {
    t_trans->clear_status();

    return t_material->clear_status();
}

int T2D2::reset_status() {
    t_trans->reset_status();

    return t_material->reset_status();
}

std::vector<vec> T2D2::record(const OutputType P) { return t_material->record(P); }

void T2D2::print() {
    suanpan_info("A 2D truss element with ");
    if(nlgeom)
        suanpan_info("corotational formulation, assuming constant {} and {} strain. ", update_area ? "volume" : "area", log_strain ? "logarithmic" : "engineering");
    else
        suanpan_info("linear formulation. ");
    suanpan_info("The nodes connected are:", node_encoding);
    suanpan_info("The area is {:.4E}. The initial element length is {:.4E}.\n", area, length);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    t_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

vtkSmartPointer<vtkCell> T2D2::GetCell() const { return vtkSmartPointer<vtkLine>::New(); }

mat T2D2::GetData(const OutputType P) {
    if(OutputType::A == P) return reshape(get_current_acceleration(), t_dof, t_node);
    if(OutputType::V == P) return reshape(get_current_velocity(), t_dof, t_node);
    if(OutputType::U == P) return reshape(get_current_displacement(), t_dof, t_node);

    vec t_stress;
    if(const auto t_data = t_material->record(P); !t_data.empty()) t_stress = t_data[0];
    return repmat(t_stress.resize(6), 1, t_node);
}

mat T2D2::GetDeformation(const double amplifier) { return get_coordinate(2).t() + amplifier * reshape(get_current_displacement(), t_dof, t_node); }

#endif
