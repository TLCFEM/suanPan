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

#include "Modifier.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

Modifier::Modifier(const unsigned T, uvec&& ET)
    : Tag(T)
    , element_tag(unique(ET)) {}

int Modifier::initialize(const shared_ptr<DomainBase>& D) {
    element_pool.clear();

    if(element_tag.empty()) {
        element_pool.reserve(D->get_element());
        for(const auto& I : D->get_element_pool()) element_pool.emplace_back(I);
    }
    else {
        element_pool.reserve(element_tag.size());
        for(const auto I : element_tag) if(D->find<Element>(I) && D->get<Element>(I)->is_active()) element_pool.emplace_back(D->get<Element>(I));
    }

    return SUANPAN_SUCCESS;
}

bool Modifier::if_apply(const shared_ptr<DomainBase>&) { return true; }

bool ModifierDynamics::if_apply(const shared_ptr<DomainBase>& D) { return AnalysisType::DYNAMICS == D->get_factory()->get_analysis_type(); }
