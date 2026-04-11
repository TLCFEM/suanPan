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
/**
 * @author tlc
 * @date 07/03/2021
 * @version 0.1.0
 * @file FixedLength.h
 * @addtogroup Constraint
 * @{
 */

#ifndef FIXEDLENGTH_H
#define FIXEDLENGTH_H

#include "Constraint.h"

#include <Domain/Factory.hpp>

/**
 * @class FixedLength
 * @brief A FixedLength class.
 *
 * The `FixedLength` constraint applies constraint to two nodes so that the
 * distance remain constant between those two nodes.
 */
template<unsigned DIM> class FixedLength : public Constraint {
    vec initial_chord;

    [[nodiscard]] bool validate_node() const final { return true; }
    [[nodiscard]] bool collect_node() const final { return true; }

protected:
    bool min_bound = false, max_bound = false;
    double min_gap = 0., max_gap = 0.;

public:
    FixedLength(const unsigned T, uvec&& N)
        : Constraint(T, 0, suanpan::translational(DIM), {}, 1) { target_node = std::move(N); }

    int initialize(const shared_ptr<DomainBase>& D) override {
        if(SUANPAN_SUCCESS != Constraint::initialize(D)) return SUANPAN_FAIL;

        initial_chord = D->get<Node>(target_node(1))->initial_position(DIM) - D->get<Node>(target_node(0))->initial_position(DIM);

        set_multiplier_size(0u);

        return SUANPAN_SUCCESS;
    }

    int process(const shared_ptr<DomainBase>& D) override {
        auto& W = D->get_factory();

        const uvec dof_i = target_node_dof.head(DIM), dof_j = target_node_dof.tail(DIM);

        const vec t_disp = W->get_trial_displacement()(dof_j) - W->get_trial_displacement()(dof_i);
        const vec t_chord = initial_chord + t_disp;

        if(const auto t_gap = dot(t_chord, t_chord); min_bound && max_bound) {
            if(0u == lagrangian_size && t_gap > min_gap && t_gap < max_gap) return SUANPAN_SUCCESS;

            auxiliary_load = (2. * std::sqrt(t_gap) < std::sqrt(min_gap) + std::sqrt(max_gap) ? min_gap : max_gap) - dot(initial_chord, initial_chord);
        }
        else if(min_bound) {
            if(0u == lagrangian_size && t_gap > min_gap) return SUANPAN_SUCCESS;

            auxiliary_load = min_gap - dot(initial_chord, initial_chord);
        }
        else if(max_bound) {
            if(0u == lagrangian_size && t_gap < max_gap) return SUANPAN_SUCCESS;

            auxiliary_load = max_gap - dot(initial_chord, initial_chord);
        }
        // no need as empty quantities are skipped
        // else auxiliary_load = 0.;

        set_multiplier_size(1u);

        auxiliary_stiffness.zeros(W->get_size(), lagrangian_size);
        auxiliary_resistance = 0.;
        for(auto I = 0u; I < DIM; ++I) {
            auxiliary_stiffness(dof_i(I)) = -(auxiliary_stiffness(dof_j(I)) = 2. * t_chord(I));
            auxiliary_resistance += t_disp(I) * (2. * initial_chord(I) + t_disp(I));
        }

        stiffness.zeros(target_node_dof.n_elem, target_node_dof.n_elem);
        const auto t_factor = 2. * trial_lambda(0);
        for(auto I = 0u; I < DIM; ++I) stiffness(I + DIM, I) = stiffness(I, I + DIM) = -(stiffness(I, I) = stiffness(I + DIM, I + DIM) = t_factor);

        resistance = auxiliary_stiffness * trial_lambda;

        return SUANPAN_SUCCESS;
    }

    [[nodiscard]] bool is_connected() const override { return true; }

    void update_status(const vec& incre_lambda) override { trial_lambda += incre_lambda; }

    void commit_status() override {
        current_lambda = trial_lambda;
        set_multiplier_size(0u);
    }

    void clear_status() override {
        current_lambda = trial_lambda.zeros();
        set_multiplier_size(0u);
    }

    void reset_status() override {
        trial_lambda = current_lambda;
        set_multiplier_size(0u);
    }
};

/**
 * @class MinimumGap
 * @brief A MinimumGap class.
 */
template<unsigned DIM> class MinimumGap final : public FixedLength<DIM> {
public:
    MinimumGap(const unsigned T, const double M, uvec&& N)
        : FixedLength<DIM>(T, std::move(N)) {
        FixedLength<DIM>::min_bound = true;
        FixedLength<DIM>::min_gap = M * M;
    }
};

/**
 * @class MaximumGap
 * @brief A MaximumGap class.
 */
template<unsigned DIM> class MaximumGap final : public FixedLength<DIM> {
public:
    MaximumGap(const unsigned T, const double M, uvec&& N)
        : FixedLength<DIM>(T, std::move(N)) {
        FixedLength<DIM>::max_bound = true;
        FixedLength<DIM>::max_gap = M * M;
    }
};

/**
 * @class Sleeve
 * @brief A Sleeve class.
 */
template<unsigned DIM> class Sleeve final : public FixedLength<DIM> {
public:
    Sleeve(const unsigned T, const double M1, const double M2, uvec&& N)
        : FixedLength<DIM>(T, std::move(N)) {
        FixedLength<DIM>::min_bound = true;
        FixedLength<DIM>::max_bound = true;
        FixedLength<DIM>::min_gap = std::pow(std::min(M1, M2), 2.);
        FixedLength<DIM>::max_gap = std::pow(std::max(M1, M2), 2.);
    }
};

/**
 * @class MaxForce
 * @brief A MaxForce class.
 */
template<unsigned DIM> class MaxForce final : public FixedLength<DIM> {
    const double max_force;

    bool trial_flag = false, current_flag = false;

public:
    MaxForce(const unsigned T, const double MF, uvec&& N)
        : FixedLength<DIM>(T, std::move(N))
        , max_force(MF) {}

    int process(const shared_ptr<DomainBase>& D) override {
        if(current_flag) {
            // if already exceeded, the constraint is not triggered
            FixedLength<DIM>::set_multiplier_size(0u);
            return SUANPAN_SUCCESS;
        }

        if(SUANPAN_SUCCESS != FixedLength<DIM>::process(D)) return SUANPAN_FAIL;

        if(0u == FixedLength<DIM>::lagrangian_size) return SUANPAN_SUCCESS;

        vec nodal_resistance(DIM);
        for(auto I = 0u; I < DIM; ++I) nodal_resistance(I) = FixedLength<DIM>::resistance(FixedLength<DIM>::target_node_dof(I));

        if(norm(nodal_resistance) > max_force) {
            trial_flag = true;
            FixedLength<DIM>::set_multiplier_size(0u);
        }

        return SUANPAN_SUCCESS;
    }

    void commit_status() override {
        current_flag = trial_flag;
        FixedLength<DIM>::commit_status();
    }

    void clear_status() override {
        current_flag = trial_flag = false;
        FixedLength<DIM>::clear_status();
    }

    void reset_status() override {
        trial_flag = current_flag;
        FixedLength<DIM>::reset_status();
    }
};

#endif

//! @}
