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

#include "CP7.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>

CP7::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& PNPXY)
    : coor(std::forward<vec>(C))
    , weight(W)
    , m_material(std::forward<unique_ptr<Material>>(M))
    , pn_pxy(std::forward<mat>(PNPXY))
    , strain_mat(3, m_size, fill::zeros) {}

CP7::CP7(const unsigned T, uvec&& N, const unsigned M, const double TH, const bool F)
    : MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, F, {DOF::U1, DOF::U2})
    , thickness(TH) {}

int CP7::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == static_cast<PlaneType>(material_proto->get_parameter(ParameterType::PLANETYPE))) suanpan::hacker(thickness) = 1.;

    const auto ele_coor = get_coordinate(m_dof);

    access::rw(characteristic_length) = sqrt(area::shoelace(mat(ele_coor.head_rows(4))));

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const IntegrationPlan plan(2, 2, IntegrationType::IRONS);

    initial_stiffness.zeros(m_size, m_size);

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(unsigned I = 0; I < plan.n_rows; ++I) {
        vec t_vec{plan(I, 0), plan(I, 1)};
        const auto pn = compute_shape_function(t_vec, 1);
        const mat jacob = pn * ele_coor;
        int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob), material_proto->get_copy(), solve(jacob, pn));

        auto& c_int_pt = int_pt.back();

        for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
            c_int_pt.strain_mat(0, K) = c_int_pt.strain_mat(2, L) = c_int_pt.pn_pxy(0, J);
            c_int_pt.strain_mat(2, K) = c_int_pt.strain_mat(1, L) = c_int_pt.pn_pxy(1, J);
        }
        initial_stiffness += c_int_pt.weight * thickness * c_int_pt.strain_mat.t() * ini_stiffness * c_int_pt.strain_mat;
    }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(const auto t_density = material_proto->get_parameter(ParameterType::DENSITY); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        for(const auto& I : int_pt) {
            const auto n_int = compute_shape_function(I.coor, 0);
            const auto t_factor = t_density * I.weight * thickness;
            for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
        }
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    body_force.zeros(m_size, m_dof);
    for(const auto& I : int_pt) {
        const mat n_int = I.weight * thickness * compute_shape_function(I.coor, 0);
        for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = 0llu; K < m_dof; ++K) body_force(L + K, K) += n_int(J);
    }

    return SUANPAN_SUCCESS;
}

int CP7::update_status() {
    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    if(nlgeom) {
        trial_geometry.zeros(m_size, m_size);

        const mat ele_disp = reshape(get_trial_displacement(), m_dof, m_node);

        mat BN(3, m_size);
        for(const auto& I : int_pt) {
            const mat gradient = ele_disp * I.pn_pxy.t() + eye(m_dof, m_dof);
            for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
                BN(0, K) = I.pn_pxy(0, J) * gradient(0, 0);
                BN(1, K) = I.pn_pxy(1, J) * gradient(0, 1);
                BN(0, L) = I.pn_pxy(0, J) * gradient(1, 0);
                BN(1, L) = I.pn_pxy(1, J) * gradient(1, 1);
                BN(2, K) = I.pn_pxy(0, J) * gradient(0, 1) + I.pn_pxy(1, J) * gradient(0, 0);
                BN(2, L) = I.pn_pxy(0, J) * gradient(1, 1) + I.pn_pxy(1, J) * gradient(1, 0);
            }

            if(I.m_material->update_trial_status(tensor::strain::to_voigt(tensor::strain::to_green(gradient))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

            const auto t_weight = I.weight * thickness;

            auto& t_stress = I.m_material->get_trial_stress();

            const auto sigma = tensor::stress::to_tensor(t_stress);

            for(unsigned J = 0, L = 0, M = 1; J < m_node; ++J, L += m_dof, M += m_dof) {
                const vec t_vec = sigma * I.pn_pxy.col(J);
                auto t_factor = t_weight * dot(I.pn_pxy.col(J), t_vec);
                trial_geometry(L, L) += t_factor;
                trial_geometry(M, M) += t_factor;
                for(auto K = J + 1, N = L + m_dof, P = N + 1; K < m_node; ++K, N += m_dof, P += m_dof) {
                    t_factor = t_weight * dot(I.pn_pxy.col(K), t_vec);
                    trial_geometry(N, L) += t_factor;
                    trial_geometry(L, N) += t_factor;
                    trial_geometry(P, M) += t_factor;
                    trial_geometry(M, P) += t_factor;
                }
            }

            trial_stiffness += t_weight * BN.t() * I.m_material->get_trial_stiffness() * BN;
            trial_resistance += t_weight * BN.t() * t_stress;
        }
    }
    else
        for(const auto& I : int_pt) {
            vec t_strain(3, fill::zeros);
            for(unsigned J = 0; J < m_node; ++J) {
                const auto& t_disp = node_ptr[J].lock()->get_trial_displacement();
                t_strain(0) += t_disp(0) * I.pn_pxy(0, J);
                t_strain(1) += t_disp(1) * I.pn_pxy(1, J);
                t_strain(2) += t_disp(0) * I.pn_pxy(1, J) + t_disp(1) * I.pn_pxy(0, J);
            }
            if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

            trial_stiffness += I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat * I.weight * thickness;
            trial_resistance += I.strain_mat.t() * I.m_material->get_trial_stress() * I.weight * thickness;
        }

    return SUANPAN_SUCCESS;
}

int CP7::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int CP7::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int CP7::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

mat CP7::compute_shape_function(const mat& coordinate, const unsigned order) const { return shape::quad(coordinate, order, m_node); }

vector<vec> CP7::record(const OutputType P) {
    vector<vec> output;
    for(const auto& I : int_pt) for(const auto& J : I.m_material->record(P)) output.emplace_back(J);
    return output;
}

void CP7::print() {
    suanpan_info("A CP7 element{}.\n", nlgeom ? " with nonlinear geometry on" : "");
    suanpan_info("The nodes connected are:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(const auto& I : int_pt) I.m_material->print();
}
