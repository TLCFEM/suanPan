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

#include "ParticleCollision.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/NodeHelper.hpp>

vec ParticleCollision::get_position(const shared_ptr<Node>& node) const {
    if(2 == num_dof) return get_trial_position<DOF::U1, DOF::U2>(node);

    return get_trial_position<DOF::U1, DOF::U2, DOF::U3>(node);
}

ParticleCollision::ParticleCollision(const unsigned T, const unsigned S, const unsigned D)
    : Constraint(T, S, 0, {}, {}, 0)
    , num_dof(D) {}

int ParticleCollision::initialize(const shared_ptr<DomainBase>& D) {
    if(const auto t_scheme = D->get_factory()->get_storage_scheme(); StorageScheme::FULL != t_scheme && StorageScheme::SPARSE != t_scheme && StorageScheme::SPARSESYMM != t_scheme) {
        suanpan_warning("The full or sparse matrix storage scheme is required.\n");
        return SUANPAN_FAIL;
    }

    return Constraint::initialize(D);
}

void ParticleCollision::apply_contact(const shared_ptr<DomainBase>& D, const shared_ptr<Node>& node_i, const shared_ptr<Node>& node_j, const bool full) {
    const auto& dof_i = node_i->get_reordered_dof();
    const auto& dof_j = node_j->get_reordered_dof();

    vec diff_pos = get_position(node_j) - get_position(node_i);
    const auto diff_norm = norm(diff_pos);

    diff_pos /= diff_norm;

    const auto force = compute_f(diff_norm);
    {
        std::scoped_lock resistance_lock(resistance_mutex);
        for(auto I = 0llu; I < diff_pos.n_elem; ++I) {
            resistance(dof_i(I)) += force * diff_pos(I);
            resistance(dof_j(I)) -= force * diff_pos(I);
        }
    }

    if(!full) return;

    auto& W = D->get_factory();
    auto& t_stiff = W->get_stiffness();

    const mat d_norm = (compute_df(diff_norm) - force / diff_norm) * diff_pos * diff_pos.t() + force / diff_norm * eye(num_dof, num_dof);

    std::scoped_lock stiffness_lock(W->get_stiffness_mutex());
    for(auto L = 0u; L < num_dof; ++L)
        for(auto K = 0u; K < num_dof; ++K) {
            t_stiff->at(dof_i(K), dof_i(L)) -= d_norm(K, L);
            t_stiff->at(dof_j(K), dof_j(L)) -= d_norm(K, L);
            t_stiff->at(dof_i(K), dof_j(L)) += d_norm(K, L);
            t_stiff->at(dof_j(K), dof_i(L)) += d_norm(K, L);
        }
}

int ParticleCollision::process(const shared_ptr<DomainBase>& D) { return process_meta(D, true); }

int ParticleCollision::process_resistance(const shared_ptr<DomainBase>& D) { return process_meta(D, false); }
