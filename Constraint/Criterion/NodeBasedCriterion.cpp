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

#include "NodeBasedCriterion.h"

#include <Domain/DomainBase.h>

NodeBasedCriterion::NodeBasedCriterion(const unsigned T, const unsigned ST, const unsigned NT, const unsigned DT, const double MA)
    : Criterion(T, ST)
    , node(NT)
    , dof(DT)
    , limit(MA) {}

int NodeBasedCriterion::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find<Node>(node)) {
        suanpan_error("Node {} is not active.\n", node);
        D->disable_criterion(get_tag());
    }

    return SUANPAN_SUCCESS;
}
