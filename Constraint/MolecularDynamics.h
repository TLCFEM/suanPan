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
 * @class MolecularDynamics
 * @brief A MolecularDynamics class.
 *
 * @author tlc
 * @date 02/12/2025
 * @version 0.1.0
 * @file MolecularDynamics.h
 * @addtogroup Constraint
 * @{
 */

#ifndef MOLECULARDYNAMICS_H
#define MOLECULARDYNAMICS_H

#include "Constraint.h"

#include <Constraint/Interaction/Interaction.h>
#include <Domain/Factory.hpp>
#include <Element/Element.h>

class MolecularDynamics : public Constraint {
    uvec interaction_tags;
    std::vector<shared_ptr<Interaction>> interactions;

    [[nodiscard]] bool validate_node() const final { return true; }

protected:
    void apply_interaction(const shared_ptr<InteractionPair>& pair) const {
        for(auto&& interaction : interactions) interaction->apply(pair);
    }

public:
    MolecularDynamics(const unsigned T, uvec&& IT, std::vector<Node::DOF>&& DC)
        : Constraint(T, 0, std::move(DC), {}, 0)
        , interaction_tags(std::move(IT)) {}

    int initialize(const shared_ptr<DomainBase>& D) override {
        if(const auto t_scheme = D->get_factory()->get_storage_scheme(); StorageScheme::FULL != t_scheme && StorageScheme::SPARSE != t_scheme && StorageScheme::SPARSESYMM != t_scheme) {
            suanpan_warning("The full or sparse matrix storage scheme is required.\n");
            return SUANPAN_FAIL;
        }

        interactions = D->get<Interaction>(interaction_tags);

        return Constraint::initialize(D);
    }
};

class MolecularDynamics2D final : public MolecularDynamics {
    struct CellList {
        int x = 0, y = 0;
        unsigned tag = 0;

        CellList(const int in_x, const int in_y, const unsigned in_tag)
            : x(in_x)
            , y(in_y)
            , tag(in_tag) {}
    };

public:
    MolecularDynamics2D(const unsigned T, uvec&& IT)
        : MolecularDynamics(T, std::move(IT), {Node::DOF::U1, Node::DOF::U2}) {}

    int process(const shared_ptr<DomainBase>& D) override {
        auto& W = D->get_factory();

        auto& element_pool = D->get_element_pool();

        suanpan::vector<CellList> list;
        list.reserve(element_pool.size());

        const auto space = 2. * std::transform_reduce(element_pool.cbegin(), element_pool.cend(), 0., [](const double a, const double b) { return std::max(a, b); }, [](const std::shared_ptr<Element>& element) { return element->get(Element::Parameter::RADIUS); });
        suanpan::for_all(element_pool, [&](const shared_ptr<Element>& element) {
            if(norm(element->get_trial_velocity()) * W->get_incre_time() > space) suanpan_warning("The nodal speed seems to be too large.\n");
            const vec new_pos = element->get_coordinate().t() + element->get_trial_displacement();
            list.emplace_back(static_cast<int>(floor(new_pos(0) / space)), static_cast<int>(floor(new_pos(1) / space)), element->get_tag());
        });

        suanpan_sort(list.begin(), list.end(), [](const CellList& a, const CellList& b) { return a.x < b.x || a.x == b.x && a.y < b.y; });

        suanpan::for_each(list.size(), [&](const size_t I) {
            for(auto J = I + 1; J < list.size(); ++J) {
                const auto diff_x = list[J].x - list[I].x;
                if(diff_x > 1) break;
                const auto diff_y = list[J].y - list[I].y;
                if(diff_x == 1 && diff_y > 1) break;
                if(std::abs(diff_y) > 1) continue;
                apply_interaction(std::make_unique<InteractionPair>(D->get<Element>(list[I].tag), D->get<Element>(list[J].tag)));
            }
        });

        return SUANPAN_SUCCESS;
    }
};

class MolecularDynamics3D final : public MolecularDynamics {
    struct CellList {
        int x = 0, y = 0, z = 0;
        unsigned tag = 0;

        CellList(const int in_x, const int in_y, const int in_z, const unsigned in_tag)
            : x(in_x)
            , y(in_y)
            , z(in_z)
            , tag(in_tag) {}
    };

public:
    MolecularDynamics3D(const unsigned T, uvec&& IT)
        : MolecularDynamics(T, std::move(IT), {Node::DOF::U1, Node::DOF::U2, Node::DOF::U3}) {}

    int process(const shared_ptr<DomainBase>& D) override {
        auto& W = D->get_factory();

        auto& element_pool = D->get_element_pool();

        suanpan::vector<CellList> list;
        list.reserve(element_pool.size());

        const auto space = 2. * std::transform_reduce(element_pool.cbegin(), element_pool.cend(), 0., [](const double a, const double b) { return std::max(a, b); }, [](const std::shared_ptr<Element>& element) { return element->get(Element::Parameter::RADIUS); });
        suanpan::for_all(element_pool, [&](const shared_ptr<Element>& element) {
            if(norm(element->get_trial_velocity()) * W->get_incre_time() > space) suanpan_warning("The nodal speed seems to be too large.\n");
            const vec new_pos = element->get_coordinate().t() + element->get_trial_displacement();
            list.emplace_back(static_cast<int>(floor(new_pos(0) / space)), static_cast<int>(floor(new_pos(1) / space)), static_cast<int>(floor(new_pos(2) / space)), element->get_tag());
        });

        suanpan_sort(list.begin(), list.end(), [](const CellList& a, const CellList& b) { return a.x < b.x || a.x == b.x && a.y < b.y || a.x == b.x && a.y == b.y && a.z < b.z; });

        suanpan::for_each(list.size(), [&](const size_t I) {
            for(auto J = I + 1; J < list.size(); ++J) {
                const auto diff_x = list[J].x - list[I].x;
                if(diff_x > 1) break;
                const auto diff_y = list[J].y - list[I].y;
                const auto diff_z = list[J].z - list[I].z;
                if(diff_x == 1 && (diff_y > 1 || diff_z > 1)) break;
                if(std::abs(diff_y) > 1 || std::abs(diff_z) > 1) continue;
                apply_interaction(std::make_unique<InteractionPair>(D->get<Element>(list[I].tag), D->get<Element>(list[J].tag)));
            }
        });

        return SUANPAN_SUCCESS;
    }
};

#endif

//! @}
