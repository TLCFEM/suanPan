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

#include "LogicCriterion.h"

#include <Domain/DomainBase.h>

LogicCriterion::LogicCriterion(const unsigned T, std::vector<unsigned>&& TT)
    : Criterion(T)
    , tags(std::move(TT)) {}

int LogicCriterion::initialize(const shared_ptr<DomainBase>& D) {
    criteria.clear();
    criteria.reserve(tags.size());

    for(const auto t_tag : tags) {
        const auto& t_criterion = D->get<Criterion>(t_tag);
        if(!t_criterion) {
            suanpan_error("Cannot find criterion {}.\n", t_tag);
            D->disable_criterion(get_tag());
            return SUANPAN_SUCCESS;
        }

        if(const auto& t_copy = criteria.emplace_back(t_criterion->unique_copy()); SUANPAN_SUCCESS != t_copy->initialize(D)) {
            suanpan_error("Fail to initialize criterion {}.\n", t_tag);
            D->disable_criterion(get_tag());
            return SUANPAN_SUCCESS;
        }
    }

    return SUANPAN_SUCCESS;
}

int LogicCriterion::process(const shared_ptr<DomainBase>& D) {
    std::vector<int> results;
    for(const auto& t_criterion : criteria)
        if(SUANPAN_FAIL == results.emplace_back(t_criterion->process(D))) return SUANPAN_FAIL;

    return check(results);
}

int LogicCriterionAll::check(const std::vector<int>& results) const {
    return std::ranges::all_of(results, [](const int result) { return SUANPAN_EXIT == result; });
}

unique_ptr<Criterion> LogicCriterionAll::unique_copy() { return std::make_unique<LogicCriterionAll>(*this); }

int LogicCriterionAny::check(const std::vector<int>& results) const {
    return std::ranges::any_of(results, [](const int result) { return SUANPAN_EXIT == result; });
}

unique_ptr<Criterion> LogicCriterionAny::unique_copy() { return std::make_unique<LogicCriterionAny>(*this); }
