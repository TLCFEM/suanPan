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

#include "RigidWallMultiplier.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

int RigidWallMultiplier::initialize(const shared_ptr<DomainBase>& D) {
    // for dynamic analysis, velocity and acceleration shall be changed
    // multiplier method does not give proper result
    if(AnalysisType::DYNAMICS == D->get_factory()->get_analysis_type()) {
        suanpan_warning("Multiplier rigid wall constraint cannot be applied in dynamic analysis.\n");
        access::rw(use_penalty) = true;
    }

    return RigidWallPenalty::initialize(D);
}

int RigidWallMultiplier::process(const shared_ptr<DomainBase>& D) {
    if(use_penalty) return RigidWallPenalty::process(D);

    auto& W = D->get_factory();

    auxiliary_stiffness.reset();
    std::vector<double> t_resistance, t_load;

    // multiplier method
    auto counter = 0llu;
    for(const auto& I : D->get_node_pool()) {
        if(!checker_handler(I)) continue;
        const vec t_pos = trial_position_handler(I) - origin;
        if(!edge_a.empty())
            if(const auto projection = dot(t_pos, edge_a); projection > length_a || projection < 0.) continue;
        if(!edge_b.empty())
            if(const auto projection = dot(t_pos, edge_b); projection > length_b || projection < 0.) continue;
        if(const auto t_pen = dot(t_pos, outer_norm); t_pen > datum::eps) continue;
        ++counter;
        auxiliary_stiffness.resize(W->get_size(), counter);
        auto& t_dof = I->get_reordered_dof();
        for(auto J = 0llu; J < n_dim; ++J) auxiliary_stiffness(t_dof(J), counter - 1) = outer_norm(J);
        const auto t_disp = trial_displacement_handler(I);
        t_resistance.emplace_back(dot(t_disp, outer_norm));
        t_load.emplace_back(-dot(t_pos - t_disp, outer_norm));
    }

    auxiliary_resistance = t_resistance;
    auxiliary_load = t_load;

    num_size = static_cast<unsigned>(auxiliary_resistance.n_elem);

    return SUANPAN_SUCCESS;
}

RigidWallMultiplier1D::RigidWallMultiplier1D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double F)
    : RigidWallMultiplier(T, S, A, resize(O, 1, 1), resize(N, 1, 1), F, 1) { set_handler<DOF::U1>(); }

RigidWallMultiplier2D::RigidWallMultiplier2D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double F)
    : RigidWallMultiplier(T, S, A, resize(O, 2, 1), resize(N, 2, 1), F, 2) { set_handler<DOF::U1, DOF::U2>(); }

RigidWallMultiplier2D::RigidWallMultiplier2D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double F)
    : RigidWallMultiplier(T, S, A, resize(O, 2, 1), resize(E1, 3, 1), resize(E2, 3, 1), F, 2) {
    set_handler<DOF::U1, DOF::U2>();
    access::rw(outer_norm).resize(2);
    access::rw(edge_a).resize(2);
    access::rw(edge_b).reset();
}

RigidWallMultiplier3D::RigidWallMultiplier3D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& N, const double F)
    : RigidWallMultiplier(T, S, A, resize(O, 3, 1), resize(N, 3, 1), F, 3) { set_handler<DOF::U1, DOF::U2, DOF::U3>(); }

RigidWallMultiplier3D::RigidWallMultiplier3D(const unsigned T, const unsigned S, const unsigned A, vec&& O, vec&& E1, vec&& E2, const double F)
    : RigidWallMultiplier(T, S, A, resize(O, 3, 1), resize(E1, 3, 1), resize(E2, 3, 1), F, 3) { set_handler<DOF::U1, DOF::U2, DOF::U3>(); }
