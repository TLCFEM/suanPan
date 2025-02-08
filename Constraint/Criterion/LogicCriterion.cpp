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

#include "LogicCriterion.h"
#include <Domain/DomainBase.h>

LogicCriterion::LogicCriterion(const unsigned T, const unsigned ST, const unsigned TA, const unsigned TB)
    : Criterion(T, ST)
    , tag_a(TA)
    , tag_b(TB) {}

int LogicCriterion::initialize(const shared_ptr<DomainBase>& D) {
    const auto& t_criterion_a = D->get<Criterion>(tag_a);
    const auto& t_criterion_b = D->get<Criterion>(tag_b);

    if(nullptr == t_criterion_a || nullptr == t_criterion_b) {
        suanpan_error("Cannot find criteria {} and/or {}.\n", tag_a, tag_b);
        D->disable_criterion(get_tag());
        return SUANPAN_SUCCESS;
    }

    criterion_a = t_criterion_a->get_copy();
    criterion_b = t_criterion_b->get_copy();

    if(SUANPAN_SUCCESS != criterion_a->initialize(D) || SUANPAN_SUCCESS != criterion_b->initialize(D)) {
        suanpan_error("Fail to initialize criteria {} and/or {}.\n", tag_a, tag_b);
        D->disable_criterion(get_tag());
    }

    return SUANPAN_SUCCESS;
}

unique_ptr<Criterion> LogicCriterionAND::get_copy() { return make_unique<LogicCriterionAND>(*this); }

int LogicCriterionAND::process(const shared_ptr<DomainBase>& D) {
    const auto result_a = criterion_a->process(D);
    const auto result_b = criterion_b->process(D);

    if(SUANPAN_FAIL == result_a || SUANPAN_FAIL == result_b) return SUANPAN_FAIL;

    return SUANPAN_EXIT == result_a && SUANPAN_EXIT == result_b;
}

unique_ptr<Criterion> LogicCriterionOR::get_copy() { return make_unique<LogicCriterionOR>(*this); }

int LogicCriterionOR::process(const shared_ptr<DomainBase>& D) {
    const auto result_a = criterion_a->process(D);
    const auto result_b = criterion_b->process(D);

    if(SUANPAN_FAIL == result_a || SUANPAN_FAIL == result_b) return SUANPAN_FAIL;

    return SUANPAN_EXIT == result_a || SUANPAN_EXIT == result_b;
}
