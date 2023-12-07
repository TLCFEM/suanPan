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

#include "PatchQuad.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/IGA/NURBSSurface.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/utility.h>

PatchQuad::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , weight(W)
    , m_material(std::move(M)) {}

PatchQuad::PatchQuad(const unsigned T, vec&& KX, vec&& KY, uvec&& N, const unsigned M, const double TH)
    : MaterialPatch2D(T, m_dof, std::move(N), uvec{M}, {KX, KY}, false)
    , m_node(static_cast<unsigned>(node_encoding.n_elem))
    , m_size(m_node * m_dof)
    , thickness(TH) {}

int PatchQuad::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    if(PlaneType::E == material_proto->get_plane_type()) suanpan::hacker(thickness) = 1.;

    const NURBSSurface2D net(knot_pool[0], knot_pool[1]);

    const auto net_size = net.get_number_of_control_points();

    field<vec> polygon(net_size(0), net_size(1));

    const auto ele_coor = get_coordinate(3);

    for(auto I = 0llu; I < prod(net_size); ++I) polygon(I) = ele_coor.row(I).t();

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const auto t_density = material_proto->get_density();

    body_force.zeros(m_size, m_dof);
    initial_stiffness.zeros(m_size, m_size);
    if(t_density > 0.) initial_mass.zeros(m_size, m_size);

    const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

    const auto ele_span = net.get_all_element_span();

    auto num_element = 1llu;
    for(auto& I : ele_span) num_element *= I.n_elem;

    int_pt.clear();
    int_pt.reserve(num_element * plan.n_rows);
    for(auto I : ele_span(0))
        for(auto J : ele_span(1)) {
            const auto &xl = knot_pool[0](I), &xh = knot_pool[0](I + 1);
            const auto &yl = knot_pool[1](J), &yh = knot_pool[1](J + 1);
            const auto dx = xh - xl, dy = yh - yl;
            for(auto K = 0u; K < plan.n_rows; ++K) {
                const auto x = xl + dx * (.5 * plan(K, 0) + .5);
                const auto y = yl + dy * (.5 * plan(K, 1) + .5);
                const auto ders = net.evaluate_shape_function_derivative(x, y, polygon, 1, 1);
                const mat pn = join_cols(.5 * dx * vectorise(ders(1, 0)).t(), .5 * dy * vectorise(ders(0, 1)).t());
                const mat jacob = pn * ele_coor.head_cols(2);
                int_pt.emplace_back(vec{x, y}, plan(K, 2) * det(jacob), material_proto->get_copy());

                auto& c_pt = int_pt.back();

                const mat pn_pxy = solve(jacob, pn);

                auto t_factor = c_pt.weight * thickness;

                c_pt.strain_mat.zeros(3, m_size);
                for(unsigned L = 0, M = 0, N = 1; L < m_node; ++L, M += m_dof, N += m_dof) {
                    c_pt.strain_mat(0, M) = c_pt.strain_mat(2, N) = pn_pxy(0, L);
                    c_pt.strain_mat(2, M) = c_pt.strain_mat(1, N) = pn_pxy(1, L);
                }
                initial_stiffness += t_factor * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;

                for(auto M = 0u, L = 0u; M < m_node; ++M, L += m_dof) for(auto N = 0llu; N < m_dof; ++N) body_force(L + N, N) += t_factor * ders(0, 0)(M);

                if(t_density > 0.) {
                    t_factor *= t_density;
                    for(auto M = 0llu; M < m_node; ++M) for(auto N = M; N < m_node; ++N) initial_mass(m_dof * M, m_dof * N) += t_factor * ders(0, 0)(M) * ders(0, 0)(N);
                }
            }
        }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(t_density > 0.) {
        for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
            initial_mass(K, K) = initial_mass(I, I);
            for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    return SUANPAN_SUCCESS;
}

int PatchQuad::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness.zeros(m_size, m_size);
    trial_resistance.zeros(m_size);

    for(const auto& I : int_pt) {
        if(I.m_material->update_trial_status(I.strain_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness += thickness * I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += thickness * I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int PatchQuad::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->commit_status();
    return code;
}

int PatchQuad::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->clear_status();
    return code;
}

int PatchQuad::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.m_material->reset_status();
    return code;
}

vector<vec> PatchQuad::record(const OutputType P) {
    vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.m_material->record(P));
    return data;
}

void PatchQuad::print() {
    suanpan_info("A four-node membrane patch (PatchQuad).\n");
    suanpan_info("The nodes connected are:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(size_t I = 0; I < int_pt.size(); ++I) {
        suanpan_info("IP {}:\t", I + 1);
        suanpan_info(int_pt[I].coor);
        int_pt[I].m_material->print();
    }
}
