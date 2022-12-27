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

#include "C3D8I.h"
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>

C3D8I::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& P)
    : coor(std::forward<vec>(C))
    , weight(W)
    , c_material(std::forward<unique_ptr<Material>>(M))
    , pn_pxyz(std::forward<mat>(P))
    , B1(6, c_size, fill::zeros)
    , B2(6, 9, fill::zeros) {}

C3D8I::C3D8I(const unsigned T, uvec&& N, const unsigned M)
    : MaterialElement3D(T, c_node, c_dof, std::forward<uvec>(N), uvec{M}, false) {}

int C3D8I::initialize(const shared_ptr<DomainBase>& D) {
    auto& mat_proto = D->get<Material>(material_tag(0));

    const auto ele_coor = get_coordinate(c_dof);

    auto& mat_stiff = mat_proto->get_initial_stiffness();

    const IntegrationPlan plan(3, 2, IntegrationType::IRONS);

    initial_stiffness.zeros(c_size, c_size);
    mat stiff_a(9, 9, fill::zeros), stiff_b(c_size, 9, fill::zeros);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1), plan(I, 2)};
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(std::move(t_vec), plan(I, 3) * det(jacob), mat_proto->get_copy(), solve(jacob, pn));

        auto& c_pt = int_pt.back();
        for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += c_dof, L += c_dof, M += c_dof) {
            c_pt.B1(0, K) = c_pt.B1(3, L) = c_pt.B1(5, M) = c_pt.pn_pxyz(0, J);
            c_pt.B1(3, K) = c_pt.B1(1, L) = c_pt.B1(4, M) = c_pt.pn_pxyz(1, J);
            c_pt.B1(5, K) = c_pt.B1(4, L) = c_pt.B1(2, M) = c_pt.pn_pxyz(2, J);
        }

        const vec pbn_pxyz = solve(jacob, -2. * c_pt.coor);
        c_pt.B2(0, 0) = c_pt.B2(3, 1) = c_pt.B2(5, 2) = pbn_pxyz(0);
        c_pt.B2(1, 4) = c_pt.B2(3, 3) = c_pt.B2(4, 5) = pbn_pxyz(1);
        c_pt.B2(2, 8) = c_pt.B2(4, 7) = c_pt.B2(5, 6) = pbn_pxyz(2);

        initial_stiffness += c_pt.weight * c_pt.B1.t() * mat_stiff * c_pt.B1;

        const auto t_stiff = mat_stiff * c_pt.B2 * c_pt.weight;
        stiff_a += c_pt.B2.t() * t_stiff;
        stiff_b += c_pt.B1.t() * t_stiff;
    }
    initial_stiffness -= stiff_b * solve(stiff_a, stiff_b.t());
    trial_stiffness = current_stiffness = initial_stiffness;

    if(const auto t_density = mat_proto->get_parameter(ParameterType::DENSITY); t_density > 0.) {
        initial_mass.zeros(c_size, c_size);
        for(const auto& I : int_pt) {
            const auto n_int = compute_shape_function(I.coor, 0);
            const auto t_factor = t_density * I.weight;
            for(auto J = 0u, L = 0u; J < c_node; ++J, L += c_dof) for(auto K = J, M = L; K < c_node; ++K, M += c_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
        }
        for(unsigned I = 0, K = 1, L = 2; I < c_size; I += c_dof, K += c_dof, L += c_dof) {
            initial_mass(K, K) = initial_mass(L, L) = initial_mass(I, I);
            for(auto J = I + c_dof, M = J + 1, N = J + 2; J < c_size; J += c_dof, M += c_dof, N += c_dof) initial_mass(J, I) = initial_mass(K, M) = initial_mass(L, N) = initial_mass(M, K) = initial_mass(N, L) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(c_size, c_dof);
    for(const auto& I : int_pt) {
        const mat n_int = I.weight * compute_shape_function(I.coor, 0);
        for(auto J = 0u, L = 0u; J < c_node; ++J, L += c_dof) for(auto K = 0llu; K < c_dof; ++K) body_force(L + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int C3D8I::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness.zeros(c_size, c_size);
    trial_resistance.zeros(c_size);
    mat stiff_a(9, 9, fill::zeros), stiff_b(c_size, 9, fill::zeros);
    vec resistance_a(9, fill::zeros);

    for(const auto& I : int_pt) {
        if(I.c_material->update_trial_status(I.B1 * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        trial_stiffness += I.weight * I.B1.t() * I.c_material->get_trial_stiffness() * I.B1;
        trial_resistance += I.weight * I.B1.t() * I.c_material->get_trial_stress();
        const auto t_stiff = I.c_material->get_trial_stiffness() * I.B2 * I.weight;
        stiff_a += I.B2.t() * t_stiff;
        stiff_b += I.B1.t() * t_stiff;
        resistance_a += I.weight * I.B2.t() * I.c_material->get_trial_stress();
    }
    trial_stiffness -= stiff_b * solve(stiff_a, stiff_b.t());
    trial_resistance -= stiff_b * solve(stiff_a, resistance_a);

    return SUANPAN_SUCCESS;
}

int C3D8I::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->commit_status();
    return code;
}

int C3D8I::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->clear_status();
    return code;
}

int C3D8I::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->reset_status();
    return code;
}

mat C3D8I::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::cube(coordinate, order, c_node); }

vector<vec> C3D8I::record(const OutputType T) {
    vector<vec> data;
    for(const auto& I : int_pt) for(const auto& J : I.c_material->record(T)) data.emplace_back(J);
    return data;
}

void C3D8I::print() {
    node_encoding.t().print("C3D8I element connects nodes:");
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(const auto& t_pt : int_pt) {
        t_pt.c_material->print();
        suanpan_info("Strain:\t");
        t_pt.c_material->get_trial_strain().t().print();
        suanpan_info("Stress:\t");
        t_pt.c_material->get_trial_stress().t().print();
    }
}

#ifdef SUANPAN_VTK
#include <vtkHexahedron.h>

void C3D8I::Setup() {
    vtk_cell = vtkSmartPointer<vtkHexahedron>::New();
    const auto ele_coor = get_coordinate(3);
    for(unsigned I = 0; I < c_node; ++I) {
        vtk_cell->GetPointIds()->SetId(I, static_cast<vtkIdType>(node_encoding(I)));
        vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
    }
}

mat C3D8I::GetData(const OutputType P) {
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

void C3D8I::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
    mat t_disp(6, c_node, fill::zeros);

    if(OutputType::A == type) t_disp.rows(0, 2) = reshape(get_current_acceleration(), c_dof, c_node);
    else if(OutputType::V == type) t_disp.rows(0, 2) = reshape(get_current_velocity(), c_dof, c_node);
    else if(OutputType::U == type) t_disp.rows(0, 2) = reshape(get_current_displacement(), c_dof, c_node);

    for(unsigned I = 0; I < c_node; ++I) arrays->SetTuple(static_cast<vtkIdType>(node_encoding(I)), t_disp.colptr(I));
}

void C3D8I::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
    const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement(), c_dof, c_node).t();
    for(unsigned I = 0; I < c_node; ++I) nodes->SetPoint(static_cast<vtkIdType>(node_encoding(I)), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
