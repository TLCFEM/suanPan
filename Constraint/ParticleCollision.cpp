/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include <Domain/Factory.hpp>

ParticleCollision::ParticleCollision(const unsigned T, const unsigned D)
    : Constraint(T, 0, suanpan::translational(D), {}, 0)
    , dimension(D) {}

int ParticleCollision::initialize(const shared_ptr<DomainBase>& D) {
    if(const auto t_scheme = D->get_factory()->get_storage_scheme(); StorageScheme::FULL != t_scheme && StorageScheme::SPARSE != t_scheme && StorageScheme::SPARSESYMM != t_scheme) {
        suanpan_warning("The full or sparse matrix storage scheme is required.\n");
        return SUANPAN_FAIL;
    }

    if(SUANPAN_SUCCESS != Constraint::initialize(D)) return SUANPAN_FAIL;
    return validate_node_impl(D) ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

void ParticleCollision::apply_contact(const shared_ptr<DomainBase>& D, const shared_ptr<Node>& node_i, const shared_ptr<Node>& node_j, const bool full) {
    const uvec dof_i = node_i->get_reordered_dof().head(dimension);
    const uvec dof_j = node_j->get_reordered_dof().head(dimension);

    vec unit_chord = node_j->trial_position(dimension) - node_i->trial_position(dimension);
    const auto distance = norm(unit_chord);

    unit_chord /= distance;

    const auto force = compute_f(distance);
    {
        std::scoped_lock resistance_lock(resistance_mutex);
        for(auto I = 0llu; I < unit_chord.n_elem; ++I) {
            resistance(dof_i(I)) += force * unit_chord(I);
            resistance(dof_j(I)) -= force * unit_chord(I);
        }
    }

    if(!full) return;

    auto& W = D->get_factory();
    auto& t_stiff = W->get_stiffness();

    const mat d_norm = (compute_df(distance) - force / distance) * unit_chord * unit_chord.t() + force / distance * eye(dimension, dimension);

    std::scoped_lock stiffness_lock(W->get_stiffness_mutex());
    for(auto L = 0u; L < dimension; ++L)
        for(auto K = 0u; K < dimension; ++K) {
            t_stiff->at(dof_i(K), dof_i(L)) -= d_norm(K, L);
            t_stiff->at(dof_j(K), dof_j(L)) -= d_norm(K, L);
            t_stiff->at(dof_i(K), dof_j(L)) += d_norm(K, L);
            t_stiff->at(dof_j(K), dof_i(L)) += d_norm(K, L);
        }
}

int ParticleCollision::process(const shared_ptr<DomainBase>& D) { return process_meta(D, true); }

int ParticleCollision::process_resistance(const shared_ptr<DomainBase>& D) { return process_meta(D, false); }
