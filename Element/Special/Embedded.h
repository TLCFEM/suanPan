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
 * @class Embedded
 * @brief A Embedded class.
 *
 * The Embedded class.
 *
 * @author tlc
 * @date 28/11/2025
 * @version 0.2.0
 * @file Embedded.h
 * @addtogroup Element
 * @{
 */

#ifndef EMBEDDED_H
#define EMBEDDED_H

#include <Domain/DomainBase.h>
#include <Element/Element.h>

template<unsigned DIM> class Embedded final : public Element {
    static constexpr unsigned max_iteration = 20u;

    const unsigned host_tag;
    const double multiplier;

    rowvec shape;
    vec weight;

    std::vector<uvec> idx;

    shared_ptr<Element> host_element;

public:
    Embedded(const unsigned T, const unsigned ET, const unsigned NT, const double P)
        : Element(T, DIM, ET, NT, translational(DIM))
        , host_tag(ET)
        , multiplier(P) {}

    int initialize(const shared_ptr<DomainBase>& D) override {
        host_element = D->get<Element>(host_tag);

        vec normalised_coor(DIM, fill::zeros);

        if(!host_element->is_active() || host_element->compute_shape_function(normalised_coor, 0).empty()) return SUANPAN_FAIL;

        const auto host_size = host_element->get_node_number();

        idx.clear();
        idx.reserve(DIM);
        idx.emplace_back(linspace<uvec>(DIM, DIM * host_size, host_size));
        for(auto I = 1u; I < DIM; ++I) idx.emplace_back(idx.back() + 1);

        const auto temp_coor = get_coordinate(DIM);
        const mat element_coor = temp_coor.tail_rows(host_size);
        const vec node_coor = temp_coor.row(0).t();

        auto counter = 0u;
        while(++counter <= max_iteration) {
            const vec incre = solve((host_element->compute_shape_function(normalised_coor, 1) * element_coor).t(), node_coor - ((shape = host_element->compute_shape_function(normalised_coor, 0)) * element_coor).t());
            if(suanpan::inf_norm(incre) < 1E-14) {
                weight.ones(get_node_number());
                weight.tail(host_size) = -shape.t();

                const vec shape_a = multiplier * shape.t();
                const mat shape_b = -shape_a * shape;

                initial_stiffness.zeros(get_total_number(), get_total_number());
                for(auto I = 0u; I < DIM; ++I) {
                    initial_stiffness(I, I) = -multiplier;
                    initial_stiffness(uvec{I}, idx[I]) = shape_a.t();
                    initial_stiffness(idx[I], uvec{I}) = shape_a;
                    initial_stiffness(idx[I], idx[I]) = shape_b;
                }
                ConstantStiffness(this);

                break;
            }
            normalised_coor += incre;
        }

        return max_iteration == counter ? SUANPAN_FAIL : SUANPAN_SUCCESS;
    }

    int update_status() override {
        trial_resistance.zeros(get_total_number());

        const vec reaction = multiplier * reshape(get_trial_displacement(), DIM, get_node_number()) * weight;

        trial_resistance.head(DIM) = -reaction;

        for(auto I = 0u; I < DIM; ++I) trial_resistance(idx[I]) = shape.t() * reaction(I);

        return SUANPAN_SUCCESS;
    }

    int clear_status() override { return SUANPAN_SUCCESS; }
    int commit_status() override { return SUANPAN_SUCCESS; }
    int reset_status() override { return SUANPAN_SUCCESS; }
};

#endif

//! @}
