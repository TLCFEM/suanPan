/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "RigidWallMultiplier.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

int RigidWallMultiplier::initialize(const shared_ptr<DomainBase>& D) {
    // for dynamic analysis, velocity and acceleration shall be changed
    // multiplier method does not give proper result
    if(AnalysisType::DYNAMICS == D->get_factory()->get_analysis_type()) {
        suanpan_warning("multiplier rigid wall constraint cannot be applied in dynamic analysis.\n");
        access::rw(use_penalty) = true;
    }

    return RigidWallPenalty::initialize(D);
}

int RigidWallMultiplier::process(const shared_ptr<DomainBase>& D) {
    if(use_penalty) return RigidWallPenalty::process(D);

    auto& W = D->get_factory();

    auxiliary_stiffness.reset();
    vector<double> t_resistance, t_load;

    // multiplier method
    auto counter = 0llu;
    for(const auto& I : D->get_node_pool()) {
        auto& t_dof = I->get_reordered_dof();
        auto& t_coor = I->get_coordinate();
        auto& t_disp = I->get_trial_displacement();
        const auto t_size = std::min(std::min(t_coor.n_elem, t_disp.n_elem), outer_norm.n_elem);
        vec t_a = -origin;
        vec t_b = zeros(size(t_a));
        for(auto J = 0llu; J < t_size; ++J) {
            t_a(J) += t_coor(J);
            t_b(J) = t_disp(J);
        }
        const vec t_c = t_a + t_b;
        if(!edge_a.empty() && dot(t_c, edge_a) > norm(edge_a) || !edge_b.empty() && dot(t_c, edge_b) > norm(edge_b)) continue;
        if(const auto t_pen = dot(t_c, outer_norm); t_pen > datum::eps) continue;
        ++counter;
        auxiliary_stiffness.resize(W->get_size(), counter);
        for(auto J = 0llu; J < t_size; ++J) auxiliary_stiffness(t_dof(J), counter - 1) = outer_norm(J);
        t_resistance.emplace_back(dot(t_b, outer_norm));
        t_load.emplace_back(-dot(t_a, outer_norm));
    }

    auxiliary_resistance = t_resistance;
    auxiliary_load = t_load;

    num_size = static_cast<unsigned>(auxiliary_resistance.n_elem);

    return SUANPAN_SUCCESS;
}
