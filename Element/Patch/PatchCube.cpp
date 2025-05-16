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

#include "PatchCube.h"

#include <Domain/DomainBase.h>
#include <Element/Utility/IGA/NURBSVolume.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>

PatchCube::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
    : coor(std::move(C))
    , weight(W)
    , c_material(std::move(M)) {}

PatchCube::PatchCube(const unsigned T, vec&& KX, vec&& KY, vec&& KZ, uvec&& N, const unsigned M)
    : MaterialPatch3D(T, c_dof, std::move(N), uvec{M}, {KX, KY, KZ}, false)
    , c_node(static_cast<unsigned>(node_encoding.n_elem)) {}

int PatchCube::initialize(const shared_ptr<DomainBase>& D) {
    auto& material_proto = D->get<Material>(material_tag(0));

    const NURBSVolume3D net(knot_pool[0], knot_pool[1], knot_pool[2]);

    const auto net_size = net.get_number_of_control_points();

    field<vec> polygon(net_size(0), net_size(1), net_size(2));

    const auto ele_coor = get_coordinate(4);

    for(auto I = 0llu; I < prod(net_size); ++I) polygon(I) = ele_coor.row(I).t();

    auto& ini_stiffness = material_proto->get_initial_stiffness();

    const auto t_density = material_proto->get_density();

    const IntegrationPlan plan(3, 2, IntegrationType::GAUSS);

    initial_stiffness.zeros(c_size, c_size);
    body_force.zeros(c_size, c_dof);
    if(t_density > 0.) initial_mass.zeros(c_size, c_size);

    const auto ele_span = net.get_all_element_span();

    auto num_element = 1llu;
    for(auto& I : ele_span) num_element *= I.n_elem;

    int_pt.clear();
    int_pt.reserve(plan.n_rows);
    for(auto I : ele_span(0))
        for(auto J : ele_span(1))
            for(auto K : ele_span(2)) {
                const auto &xl = knot_pool[0](I), &xh = knot_pool[0](I + 1);
                const auto &yl = knot_pool[1](J), &yh = knot_pool[1](J + 1);
                const auto &zl = knot_pool[2](K), &zh = knot_pool[2](K + 1);
                const auto dx = xh - xl, dy = yh - yl, dz = zh - zl;
                for(unsigned L = 0; L < plan.n_rows; ++L) {
                    const auto x = xl + dx * (.5 * plan(L, 0) + .5);
                    const auto y = yl + dy * (.5 * plan(L, 1) + .5);
                    const auto z = zl + dz * (.5 * plan(L, 2) + .5);
                    const auto ders = net.evaluate_shape_function_derivative(x, y, z, polygon, 1, 1, 1);
                    const auto pn = join_cols(.5 * dx * vectorise(ders(1, 0, 0)).t(), .5 * dy * vectorise(ders(0, 1, 0)).t(), .5 * dz * vectorise(ders(0, 0, 1)).t());
                    const mat jacob = pn * ele_coor.head_cols(3);
                    int_pt.emplace_back(vec{x, y, z}, plan(L, 3) * det(jacob), material_proto->get_copy());

                    auto& c_pt = int_pt.back();

                    const mat pn_pxyz = solve(jacob, pn);

                    c_pt.strain_mat.zeros(6, c_size);
                    for(unsigned P = 0, M = 0, N = 1, O = 2; P < c_node; ++P, M += c_dof, N += c_dof, O += c_dof) {
                        c_pt.strain_mat(0, M) = c_pt.strain_mat(3, N) = c_pt.strain_mat(5, O) = pn_pxyz(0, P);
                        c_pt.strain_mat(3, M) = c_pt.strain_mat(1, N) = c_pt.strain_mat(4, O) = pn_pxyz(1, P);
                        c_pt.strain_mat(5, M) = c_pt.strain_mat(4, N) = c_pt.strain_mat(2, O) = pn_pxyz(2, P);
                    }
                    initial_stiffness += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;

                    for(auto M = 0u, N = 0u; M < c_node; ++M, N += c_dof)
                        for(auto P = 0llu; P < c_dof; ++P) body_force(N + P, P) += c_pt.weight * ders(0, 0, 0)(M);

                    if(t_density > 0.)
                        for(auto M = 0llu; M < c_node; ++M)
                            for(auto N = M; N < c_node; ++N) initial_mass(c_dof * M, c_dof * N) += t_density * c_pt.weight * ders(0, 0, 0)(M) * ders(0, 0, 0)(N);
                }
            }
    trial_stiffness = current_stiffness = initial_stiffness;

    if(t_density > 0.) {
        for(auto I = 0u, K = 1u, L = 2u; I < c_size; I += c_dof, K += c_dof, L += c_dof) {
            initial_mass(K, K) = initial_mass(L, L) = initial_mass(I, I);
            for(auto J = I + c_dof, M = J + 1, N = J + 2; J < c_size; J += c_dof, M += c_dof, N += c_dof) initial_mass(J, I) = initial_mass(K, M) = initial_mass(L, N) = initial_mass(M, K) = initial_mass(N, L) = initial_mass(I, J);
        }
        ConstantMass(this);
    }

    return SUANPAN_SUCCESS;
}

int PatchCube::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_stiffness.zeros(c_size, c_size);
    trial_resistance.zeros(c_size);

    for(const auto& I : int_pt) {
        if(I.c_material->update_trial_status(I.strain_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        trial_stiffness += I.weight * I.strain_mat.t() * I.c_material->get_trial_stiffness() * I.strain_mat;
        trial_resistance += I.weight * I.strain_mat.t() * I.c_material->get_trial_stress();
    }

    return SUANPAN_SUCCESS;
}

int PatchCube::commit_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->commit_status();
    return code;
}

int PatchCube::clear_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->clear_status();
    return code;
}

int PatchCube::reset_status() {
    auto code = 0;
    for(const auto& I : int_pt) code += I.c_material->reset_status();
    return code;
}

std::vector<vec> PatchCube::record(const OutputType P) {
    std::vector<vec> data;
    for(const auto& I : int_pt) append_to(data, I.c_material->record(P));
    return data;
}

void PatchCube::print() {
    suanpan_info("PatchCube element connects nodes:", node_encoding);
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    for(const auto& t_pt : int_pt) {
        t_pt.c_material->print();
        suanpan_info("Strain:\t", t_pt.c_material->get_trial_strain());
        suanpan_info("Stress:\t", t_pt.c_material->get_trial_stress());
    }
}
