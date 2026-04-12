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

template<unsigned DIM, bool ROTATION> class MolecularDynamics : public Constraint {
    uvec interaction_tags;
    std::vector<shared_ptr<Interaction>> interactions;

protected:
    std::vector<shared_ptr<Element>> elements;
    double space = 0.;

    void apply_interaction(const bool full, const shared_ptr<Element>& element) const {
        for(auto&& interaction : interactions) interaction->apply(full, element);
    }

    void apply_interaction(const bool full, const shared_ptr<InteractionPair>& pair) const {
        pair->set_dimension(DIM);
        pair->set_inertial(ROTATION);
        for(auto&& interaction : interactions) interaction->apply(full, pair);
    }

    [[nodiscard]] virtual int process_impl(const shared_ptr<DomainBase>&, bool) = 0;

    int process_entrypoint(const shared_ptr<DomainBase>& D, const bool full) {
        suanpan::for_all(D->get_element_pool(), [&](const shared_ptr<Element>& element) { apply_interaction(full, element); });
        return process_impl(D, full);
    }

public:
    MolecularDynamics(const unsigned T, uvec&& IT)
        : Constraint(T, 0, ROTATION ? suanpan::mechanical(DIM) : suanpan::translational(DIM), {}, 0)
        , interaction_tags(std::move(IT)) {}

    int initialize(const shared_ptr<DomainBase>& D) override {
        if(const auto t_scheme = D->get_factory()->get_storage_scheme(); StorageScheme::FULL != t_scheme && StorageScheme::SPARSE != t_scheme && StorageScheme::SPARSESYMM != t_scheme) {
            suanpan_warning("The full or sparse matrix storage scheme is required.\n");
            return SUANPAN_FAIL;
        }

        interactions.clear();
        for(auto&& item : D->get<Interaction>(interaction_tags))
            if(item && item->is_active()) interactions.emplace_back(item);

        elements.clear();
        for(auto&& item : D->get_element_pool())
            if(item && item->is_active() && item->type() == Element::Type::DEM) elements.emplace_back(item);

        space = 2. * std::transform_reduce(elements.cbegin(), elements.cend(), 0., [](const double a, const double b) { return std::max(a, b); }, [](const std::shared_ptr<Element>& element) { return element->get(Element::Parameter::RADIUS); });

        if(SUANPAN_SUCCESS != Constraint::initialize(D)) return SUANPAN_FAIL;

        if(!validate_node_impl(D)) return SUANPAN_FAIL;

        return SUANPAN_SUCCESS;
    }

    int process(const shared_ptr<DomainBase>& D) override { return process_entrypoint(D, true); }

    int process_resistance(const shared_ptr<DomainBase>& D) override { return process_entrypoint(D, false); }
};

class MolecularDynamics2D final : public MolecularDynamics<2u, false> {
    struct CellList {
        int x = 0, y = 0;
        unsigned tag = 0;

        CellList(const int in_x, const int in_y, const unsigned in_tag)
            : x(in_x)
            , y(in_y)
            , tag(in_tag) {}
    };

    int process_impl(const shared_ptr<DomainBase>& D, const bool full) override {
        auto& W = D->get_factory();

        suanpan::vector<CellList> list;
        list.reserve(elements.size());

        suanpan::for_all(elements, [&](const shared_ptr<Element>& element) {
            if(norm(element->get_trial_velocity().head(2)) * W->get_incre_time() > space) suanpan_warning("The nodal speed seems to be too large.\n");
            const vec new_pos = element->get_coordinate(2).t() + element->get_trial_displacement().head(2);
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
                apply_interaction(full, std::make_unique<InteractionPair>(D->get<Element>(list[I].tag), D->get<Element>(list[J].tag)));
            }
        });

        return SUANPAN_SUCCESS;
    }

public:
    using MolecularDynamics ::MolecularDynamics;
};

class MolecularDynamics3D final : public MolecularDynamics<3u, false> {
    struct CellList {
        int x = 0, y = 0, z = 0;
        unsigned tag = 0;

        CellList(const int in_x, const int in_y, const int in_z, const unsigned in_tag)
            : x(in_x)
            , y(in_y)
            , z(in_z)
            , tag(in_tag) {}
    };

    int process_impl(const shared_ptr<DomainBase>& D, const bool full) override {
        auto& W = D->get_factory();

        suanpan::vector<CellList> list;
        list.reserve(elements.size());

        suanpan::for_all(elements, [&](const shared_ptr<Element>& element) {
            if(norm(element->get_trial_velocity().head(3)) * W->get_incre_time() > space) suanpan_warning("The nodal speed seems to be too large.\n");
            const vec new_pos = element->get_coordinate(3).t() + element->get_trial_displacement().head(3);
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
                apply_interaction(full, std::make_unique<InteractionPair>(D->get<Element>(list[I].tag), D->get<Element>(list[J].tag)));
            }
        });

        return SUANPAN_SUCCESS;
    }

public:
    using MolecularDynamics ::MolecularDynamics;
};

#endif

//! @}
