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

#include "LeeNewmark.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

uword LeeNewmark::get_total_size() const { return n_damping * n_block + n_block; }

void LeeNewmark::update_stiffness() const {
    if(factory->is_sparse())
        for(uword I = 0, J = n_block; I < n_damping; ++I, J += n_block) {
            stiffness->triplet_mat.assemble(current_mass->triplet_mat, {J, 0, J}, {J, J, 0}, {mass_coef(I), -mass_coef(I), -C1 * mass_coef(I)});
            stiffness->triplet_mat.assemble(current_stiffness->triplet_mat, J, J, stiffness_coef(I));
        }
    else
        for(uword I = 0, J = n_block; I < n_damping; ++I, J += n_block)
            for(unsigned K = 0; K < n_block; ++K) {
                const auto M = K + J;
                for(unsigned L = 0; L < n_block; ++L) {
                    const auto N = L + J;
                    sp_d auto t_val = current_mass->operator()(K, L);
                    if(!suanpan::approx_equal(0., t_val)) stiffness->at(M, L) = C1 * (stiffness->at(K, N) = -(stiffness->at(M, N) = mass_coef(I) * t_val));
                    t_val = current_stiffness->operator()(K, L);
                    if(!suanpan::approx_equal(0., t_val)) stiffness->at(M, N) = stiffness_coef(I) * t_val;
                }
            }
}

void LeeNewmark::update_residual() const {
    const auto& t_vel = factory->get_trial_velocity();

    for(uword I = 0, J = n_block, K = J + n_block - 1; I < n_damping; ++I, J += n_block, K += n_block) {
        const vec n_internal(&trial_internal(J), n_block);
        residual.rows(J, K) = current_mass * vec(t_vel - n_internal) * mass_coef(I) - current_stiffness * n_internal * stiffness_coef(I);
    }
}

void LeeNewmark::initialize_mass(const shared_ptr<DomainBase>&) {
    // assuming mass does not change
    // otherwise swap and assemble
    current_mass = factory->get_mass()->make_copy();
    if(if_iterative) {
        const auto current_mass_max = current_mass->max() * 1E-10;
        for(uword I = 0llu; I < std::min(current_mass->n_rows, current_mass->n_cols); ++I) current_mass->at(I, I) += current_mass_max;
    }
}

void LeeNewmark::initialize_stiffness(const shared_ptr<DomainBase>& D) {
    auto& t_stiff = factory->modify_stiffness();
    auto& t_geometry = factory->modify_geometry();

    auto fa = std::async([&] {
        current_stiffness.swap(t_stiff);
        D->assemble_current_stiffness();
    });

    if(!factory->is_nlgeom()) fa.get();
    else {
        current_geometry.swap(t_geometry);
        D->assemble_current_geometry();
        current_geometry.swap(t_geometry);
        fa.get();
        t_stiff += current_geometry;
    }
}

LeeNewmark::LeeNewmark(const unsigned T, vec&& X, vec&& F, const double A, const double B)
    : LeeNewmarkBase(T, A, B)
    , mass_coef(4. * X % F)
    , stiffness_coef(4. * X / F)
    , CM(4. * dot(X, F)) {}

int LeeNewmark::initialize() {
    if(SUANPAN_SUCCESS != LeeNewmarkBase::initialize()) return SUANPAN_FAIL;

    current_mass = factory->get_mass()->make_copy();
    current_stiffness = factory->get_stiffness()->make_copy();
    if(factory->is_nlgeom()) current_geometry = factory->get_geometry()->make_copy();

    return SUANPAN_SUCCESS;
}

int LeeNewmark::process_constraint() {
    const auto& D = get_domain();

    // process constraint for the first time to obtain proper stiffness
    if(SUANPAN_SUCCESS != LeeNewmarkBase::process_constraint()) return SUANPAN_FAIL;

    // this stiffness contains geometry, mass and damping which are handled in Newmark::assemble_matrix()
    auto& t_stiff = factory->modify_stiffness();

    t_stiff->csc_condense();

    if(first_iteration) {
        stiffness->triplet_mat.init((4 * n_damping + 2) * t_stiff->n_elem);

        initialize_mass(D);
    }
    else if(SUANPAN_SUCCESS != erase_top_left_block()) return SUANPAN_FAIL;

    t_stiff += C1 * CM * current_mass;

    t_stiff->csc_condense();

    // check in tangent stiffness
    stiffness += t_stiff;

    if(first_iteration) {
        // for the first iteration of each substep
        // store current stiffness to be used in the whole substep
        // check in constant terms that does not change in the substep
        initialize_stiffness(D);

        if(SUANPAN_SUCCESS != LeeNewmarkBase::process_constraint()) return SUANPAN_FAIL;
        t_stiff->csc_condense();
        current_stiffness.swap(t_stiff);

        update_stiffness();

        first_iteration = false;
    }

    update_residual();

    return SUANPAN_SUCCESS;
}

int LeeNewmark::process_constraint_resistance() {
    update_residual();

    return LeeNewmarkBase::process_constraint_resistance();
}

void LeeNewmark::assemble_resistance() {
    const auto& D = get_domain();
    const auto& W = factory;

    D->assemble_resistance();
    D->assemble_inertial_force();
    // consider independent viscous device
    D->assemble_damping_force();

    if(nullptr != current_mass) {
        vec internal_velocity = CM * W->get_trial_velocity();
        for(uword I = 0, J = n_block; I < n_damping; ++I, J += n_block) {
            const vec n_internal(&trial_internal(J), n_block);
            internal_velocity -= mass_coef(I) * n_internal;
        }
        W->update_trial_damping_force_by(current_mass * internal_velocity);
    }

    W->set_sushi(W->get_trial_resistance() + W->get_trial_damping_force() + W->get_trial_inertial_force());
}

void LeeNewmark::print() {
    suanpan_info("A Newmark solver using Lee's damping model. doi: 10.1016/j.jsv.2020.115312\n");
    const vec X = .25 * sqrt(mass_coef % stiffness_coef);
    const vec F = sqrt(mass_coef / stiffness_coef);
    for(auto I = 0llu; I < n_damping; ++I)
        suanpan_info("\tDamping Ratio: {:.4f}\tFrequency (rad/s): {:.4f}\n", X(I), F(I));
}

void LeeElementalNewmark::initialize_mass(const shared_ptr<DomainBase>& D) {
    auto& t_mass = factory->modify_mass();

    current_mass.swap(t_mass);
    D->assemble_mass_container();
    current_mass.swap(t_mass);

    if(if_iterative) {
        const auto current_mass_max = current_mass->max() * 1E-10;
        for(uword I = 0llu; I < std::min(current_mass->n_rows, current_mass->n_cols); ++I) current_mass->at(I, I) += current_mass_max;
    }
}

void LeeElementalNewmark::initialize_stiffness(const shared_ptr<DomainBase>& D) {
    auto& t_stiff = factory->modify_stiffness();

    current_stiffness.swap(t_stiff);
    D->assemble_stiffness_container();
}
