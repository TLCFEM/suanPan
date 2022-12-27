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

#include "DC3D8.h"
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>

const uvec DC3D8::u_dof{0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 16, 17, 18, 20, 21, 22, 24, 25, 26, 28, 29, 30};
const uvec DC3D8::d_dof{3, 7, 11, 15, 19, 23, 27, 31};

DC3D8::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& N, mat&& P)
    : coor(std::forward<vec>(C))
    , weight(W)
    , c_material(std::forward<unique_ptr<Material>>(M))
    , n_mat(std::forward<mat>(N))
    , pn_mat(std::forward<mat>(P))
    , strain_mat(6, 24, fill::zeros) {}

DC3D8::DC3D8(const unsigned T, uvec&& N, const unsigned M, const double CL, const double RR)
    : MaterialElement3D(T, c_node, c_dof, std::forward<uvec>(N), uvec{M}, false, {DOF::U1, DOF::U2, DOF::U3, DOF::DMG})
    , release_rate(RR) { access::rw(characteristic_length) = CL; }

int DC3D8::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    const auto ele_coor = get_coordinate(3);

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(3, 2, IntegrationType::GAUSS);

    initial_stiffness.zeros(c_size, c_size);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1), plan(I, 2)};
        const auto pn = shape::cube(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(std::move(t_vec), plan(I, 3) * det(jacob), material_proto->get_copy(), shape::cube(t_vec, 0), solve(jacob, pn));

        auto& c_pt = int_pt.back();
        for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += 3, L += 3, M += 3) {
            c_pt.strain_mat(0, K) = c_pt.strain_mat(3, L) = c_pt.strain_mat(5, M) = c_pt.pn_mat(0, J);
            c_pt.strain_mat(3, K) = c_pt.strain_mat(1, L) = c_pt.strain_mat(4, M) = c_pt.pn_mat(1, J);
            c_pt.strain_mat(5, K) = c_pt.strain_mat(4, L) = c_pt.strain_mat(2, M) = c_pt.pn_mat(2, J);
        }
        initial_stiffness(u_dof, u_dof) += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
        initial_stiffness(d_dof, d_dof) += c_pt.weight * release_rate / characteristic_length * c_pt.n_mat.t() * c_pt.n_mat;
        initial_stiffness(d_dof, d_dof) += c_pt.weight * release_rate * characteristic_length * c_pt.pn_mat.t() * c_pt.pn_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

int DC3D8::update_status() {
    const auto t_disp = get_trial_displacement();

    const vec t_damage = t_disp(d_dof);

    trial_stiffness.zeros(c_size, c_size);
    trial_resistance.zeros(c_size);

    for(const auto& I : int_pt) {
        if(I.c_material->update_trial_status(I.strain_mat * t_disp(u_dof)) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        const auto pow_term = 1. - dot(t_damage, I.n_mat);
        const auto damage = pow(pow_term, 2.);

        trial_stiffness(u_dof, u_dof) += I.weight * damage * I.strain_mat.t() * I.c_material->get_trial_stiffness() * I.strain_mat;
        trial_stiffness(u_dof, d_dof) -= I.weight * 2. * pow_term * I.strain_mat.t() * I.c_material->get_trial_stress() * I.n_mat;
        trial_stiffness(d_dof, d_dof) += I.weight * (2. * I.maximum_energy + release_rate / characteristic_length) * I.n_mat.t() * I.n_mat;
        trial_stiffness(d_dof, d_dof) += I.weight * release_rate * characteristic_length * I.pn_mat.t() * I.pn_mat;

        trial_resistance(u_dof) += I.weight * damage * I.strain_mat.t() * I.c_material->get_trial_stress();
        trial_resistance(d_dof) -= I.weight * 2. * I.n_mat.t() * I.maximum_energy;
    }

    trial_resistance(d_dof) += trial_stiffness(d_dof, d_dof) * t_damage;

    return SUANPAN_SUCCESS;
}

int DC3D8::commit_status() {
    auto code = 0;
    for(auto& I : int_pt) {
        I.commit_status(I.c_material);
        I.maximum_energy = std::max(I.maximum_energy, I.strain_energy);
        code += I.c_material->commit_status();
    }
    return code;
}

int DC3D8::clear_status() {
    auto code = 0;
    for(auto& I : int_pt) {
        I.clear_status();
        I.maximum_energy = 0.;
        code += I.c_material->clear_status();
    }
    return code;
}

int DC3D8::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->reset_status();
    return code;
}

vector<vec> DC3D8::record(const OutputType T) {
    if(T == OutputType::DAMAGE) return {get_current_displacement()(d_dof)};

    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.c_material->record(T)) data.emplace_back(J);
    return data;
}

void DC3D8::print() {
    node_encoding.t().print("DC3D8 element connects:");
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(const auto& t_pt : int_pt) {
        t_pt.c_material->print();
        suanpan_info("Strain:\t");
        t_pt.c_material->get_current_strain().t().print();
        suanpan_info("Stress:\t");
        t_pt.c_material->get_current_stress().t().print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkHexahedron.h>

void DC3D8::Setup() {
    vtk_cell = vtkSmartPointer<vtkHexahedron>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < c_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

mat DC3D8::GetData(const OutputType P) {
    if(P == OutputType::DAMAGE) {
        mat t_damage(6, c_node, fill::zeros);
        t_damage.row(0) = get_current_displacement()(d_dof).t();
        return t_damage;
    }

    mat A(int_pt.size(), 7);
    mat B(int_pt.size(), 6, fill::zeros);

    for(size_t I = 0; I < int_pt.size(); ++I) {
        if(const auto C = int_pt[I].c_material->record(P); !C.empty()) B(I, 0, size(C[0])) = C[0];
        A.row(I) = interpolation::linear(int_pt[I].coor);
    }

    mat data(c_node, 7);

    data.row(0) = interpolation::linear(-1., -1., -1.);
    data.row(1) = interpolation::linear(1., -1., -1.);
    data.row(2) = interpolation::linear(1., 1., -1.);
    data.row(3) = interpolation::linear(-1., 1., -1.);
    data.row(4) = interpolation::linear(-1., -1., 1.);
    data.row(5) = interpolation::linear(1., -1., 1.);
    data.row(6) = interpolation::linear(1., 1., 1.);
    data.row(7) = interpolation::linear(-1., 1., 1.);

    return (data * solve(A, B)).t();
}

void DC3D8::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, c_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 2) = reshape(get_current_acceleration()(u_dof), 3, c_node);
    else if(OutputType::V == type) t_disp.rows(0, 2) = reshape(get_current_velocity()(u_dof), 3, c_node);
    else if(OutputType::U == type) t_disp.rows(0, 2) = reshape(get_current_displacement()(u_dof), 3, c_node);

    for(unsigned I = 0; I < c_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void DC3D8::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement()(u_dof), 3, c_node).t();
    for(unsigned I = 0; I < c_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
