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
/**
 * @class Embed
 * @brief A Embed class.
 *
 * @author tlc
 * @date 12/08/2020
 * @version 0.1.0
 * @file Embed.h
 * @addtogroup Constraint
 * @{
 */

#ifndef EMBED_H
#define EMBED_H

#include "Constraint.h"

#include <Domain/Factory.hpp>
#include <Element/Element.h>

template<unsigned DIM> class Embed final : public Constraint {
    static constexpr unsigned max_iteration = 20u;
    static constexpr double tolerance = 1E-14;

    [[nodiscard]] bool validate_node() const override { return true; }
    [[nodiscard]] bool validate_element() const override { return true; }

public:
    Embed(const unsigned T, const unsigned ET, const unsigned NT)
        : Constraint(T, 0, translational(DIM), {}, DIM) {
        target_node = NT;
        target_element = ET;
    }

    int initialize(const shared_ptr<DomainBase>& D) override {
        if(SUANPAN_SUCCESS != Constraint::initialize(D)) return SUANPAN_FAIL;

        auto& t_node = D->get<Node>(target_node(0));
        auto& t_element = D->get<Element>(target_element(0));

        vec normalised_coor(DIM, fill::zeros);

        if(t_element->compute_shape_function(normalised_coor, 0).is_empty()) return SUANPAN_FAIL;

        const auto element_coor = t_element->get_coordinate(DIM);
        const auto node_coor = t_node->initial_position(DIM);

        rowvec shape;

        // solve local system to obtain the shape function at the position of node
        auto counter = 0u;
        while(true) {
            if(max_iteration == ++counter) return SUANPAN_FAIL;
            const vec incre = solve((t_element->compute_shape_function(normalised_coor, 1) * element_coor).t(), node_coor - ((shape = t_element->compute_shape_function(normalised_coor, 0)) * element_coor).t());
            if(suanpan::inf_norm(incre) < tolerance) break;
            normalised_coor += incre;
        }

        auxiliary_stiffness.zeros(D->get_factory()->get_size(), DIM);

        auto& node_dof = t_node->get_reordered_dof();
        auto& element_dof = t_element->get_dof_encoding();
        for(auto K = 0u; K < DIM; ++K) {
            auxiliary_stiffness(node_dof(K), K) = -1.;
            for(uword I = 0, J = K; I < shape.n_elem; ++I, J += t_element->get_dof_number()) auxiliary_stiffness(element_dof(J), K) = shape(I);
        }

        return SUANPAN_SUCCESS;
    }

    int process(const shared_ptr<DomainBase>& D) override {
        auxiliary_resistance = auxiliary_stiffness.t() * D->get_factory()->get_trial_displacement();

        return SUANPAN_SUCCESS;
    }
};

#endif

//! @}
