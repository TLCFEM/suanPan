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

#include "F21H.h"

#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

F21H::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
    : coor(C)
    , weight(W)
    , b_section(std::move(M))
    , B(2, 3, fill::zeros) {
    B(0, 0) = 1.;
    B(1, 1) = .5 * (coor - 1.);
    B(1, 2) = .5 * (coor + 1.);
}

F21H::F21H(const unsigned T, uvec&& N, const unsigned S, const double L, const bool F)
    : SectionElement2D(T, b_node, b_dof, std::move(N), uvec{S}, F)
    , hinge_length(L > .5 ? .5 : L)
    , b_trans(F ? std::make_unique<B2DC>() : std::make_unique<B2DL>()) {}

int F21H::initialize(const shared_ptr<DomainBase>& D) {
    auto& section_proto = D->get<Section>(section_tag(0));

    b_trans->set_element_ptr(this);

    access::rw(length) = b_trans->get_length();

    elastic_section_flexibility = inv(section_proto->get_initial_stiffness());

    // perform integration of elastic region
    const IntegrationPlan plan(1, 2, IntegrationPlan::Type::GAUSS);
    // add two inner points of Radau quadrature
    const auto int_pt_num = plan.n_rows + 2;
    const auto elastic_length = 1. - 8. * hinge_length;
    // elastic part will be reused in computation
    elastic_local_flexibility.zeros(3, 3);
    // build up the elastic interior
    elastic_int_pt.clear();
    elastic_int_pt.reserve(int_pt_num);
    for(unsigned I = 0; I < int_pt_num; ++I) {
        double coor, weight;
        if(I == 0) {
            // left inner Radau point
            coor = 16. / 3. * hinge_length - 1.;
            weight = 3. * hinge_length;
        }
        else if(I == int_pt_num - 1) {
            // right inner Radau point
            coor = 1. - 16. / 3. * hinge_length;
            weight = 3. * hinge_length;
        }
        else {
            // Gauss points
            coor = plan(I - 1, 0) * elastic_length;
            weight = .5 * plan(I - 1, 1) * elastic_length;
        }
        elastic_int_pt.emplace_back(coor, weight, section_proto->get_copy());
        elastic_local_flexibility += elastic_int_pt[I].B.t() * elastic_section_flexibility * elastic_int_pt[I].B * weight * length;
    }

    // perform integration of hinge part
    initial_local_flexibility = elastic_local_flexibility;
    int_pt.clear();
    int_pt.reserve(2);
    int_pt.emplace_back(-1., hinge_length, section_proto->get_copy());
    int_pt.emplace_back(1., hinge_length, section_proto->get_copy());
    for(auto& I : int_pt) initial_local_flexibility += I.B.t() * elastic_section_flexibility * I.B * I.weight * length;

    trial_local_flexibility = current_local_flexibility = initial_local_flexibility;

    trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(inv(initial_local_flexibility));

    if(const auto linear_density = section_proto->get_linear_density(); linear_density > 0.) trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(linear_density);

    trial_local_deformation = current_local_deformation.zeros(3);
    trial_local_resistance = current_local_resistance.zeros(3);

    return SUANPAN_SUCCESS;
}

int F21H::update_status() {
    b_trans->update_status();

    vec residual_deformation = -trial_local_deformation;
    // transform global deformation to local one (remove rigid body motion)
    trial_local_deformation = b_trans->to_local_vec(get_trial_displacement());
    // initial residual be aware of how to compute it
    residual_deformation += trial_local_deformation;

    trial_local_resistance += solve(trial_local_flexibility, residual_deformation);

    residual_deformation.zeros();
    // consider elastic interior no residual generated
    trial_local_flexibility = elastic_local_flexibility;
    // computation of hinge part
    for(const auto& I : int_pt) {
        const vec target_resistance = I.B * trial_local_resistance;
        // compute unbalanced deformation use section stiffness
        vec incre_deformation;
        if(!solve(incre_deformation, I.b_section->get_trial_stiffness(), target_resistance - I.b_section->get_trial_resistance())) return SUANPAN_FAIL;
        // update status
        if(SUANPAN_SUCCESS != I.b_section->update_trial_status(I.b_section->get_trial_deformation() + incre_deformation)) return SUANPAN_FAIL;
        // collect new flexibility and deformation
        trial_local_flexibility += I.weight * length * I.B.t() * solve(I.b_section->get_trial_stiffness(), I.B);
        residual_deformation += I.weight * length * I.B.t() * solve(I.b_section->get_trial_stiffness(), target_resistance - I.b_section->get_trial_resistance());
    }

    trial_local_resistance -= solve(trial_local_flexibility, residual_deformation);

    trial_stiffness = b_trans->to_global_stiffness_mat(inv(trial_local_flexibility));
    trial_resistance = b_trans->to_global_vec(trial_local_resistance);

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(trial_local_resistance);

    return SUANPAN_SUCCESS;
}

int F21H::clear_status() {
    b_trans->clear_status();

    trial_local_flexibility = current_local_flexibility = initial_local_flexibility;
    trial_local_deformation = current_local_deformation.zeros();
    trial_local_resistance = current_local_resistance.zeros();
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->clear_status();
    return code;
}

int F21H::commit_status() {
    b_trans->commit_status();

    current_local_flexibility = trial_local_flexibility;
    current_local_deformation = trial_local_deformation;
    current_local_resistance = trial_local_resistance;
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->commit_status();
    return code;
}

int F21H::reset_status() {
    b_trans->reset_status();

    trial_local_flexibility = current_local_flexibility;
    trial_local_deformation = current_local_deformation;
    trial_local_resistance = current_local_resistance;
    auto code = 0;
    for(const auto& I : int_pt) code += I.b_section->reset_status();
    return code;
}

std::vector<vec> F21H::record(const OutputType P) const {
    if(P == OutputType::BEAME) return {current_local_deformation};
    if(P == OutputType::BEAMS) return {current_local_resistance};

    std::vector<vec> data;
    for(const auto& I : int_pt) suanpan::append_to(data, I.b_section->record(P));
    return data;
}

void F21H::print() {
    suanpan_info("A 2D force based beam element with lumped plasticity{} doi:https://doi.org/10.1016/0045-7949(95)00103-N \n", nlgeom ? " and corotational formulation." : ".");
    suanpan_info("The element connects nodes:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Section:\n");
    auto J = 1;
    for(const auto& I : int_pt) {
        suanpan_info("IP {}: ", J++);
        I.b_section->print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

vtkSmartPointer<vtkCell> F21H::GetCell() const { return vtkSmartPointer<vtkLine>::New(); }

mat F21H::GetData(const OutputType P) {
    const auto remap = [&](vec&& in) {
        mat data(6, b_node, fill::zeros);
        data.rows(uvec{0, 1, 5}) = in.reshape(b_dof, b_node);
        return data;
    };

    if(OutputType::A == P) return remap(get_current_acceleration());
    if(OutputType::V == P) return remap(get_current_velocity());
    if(OutputType::U == P) return remap(get_current_displacement());

    vec low, high;
    if(const auto t_data = int_pt.front().b_section->record(P); !t_data.empty()) low = t_data[0];
    if(const auto t_data = int_pt.back().b_section->record(P); !t_data.empty()) high = t_data[0];

    mat data(6, b_node);
    data.col(0) = low.resize(6);
    data.col(1) = high.resize(6);

    return data;
}

mat F21H::GetDeformation(const double amplifier) { return get_coordinate(2).t() + amplifier * get_current_displacement().reshape(b_dof, b_node).head_rows(2); }

#endif
