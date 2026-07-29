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

#include "MinMaxHistory.h"

#include <Domain/DomainBase.h>
#include <Element/Element.h>

HistoryCriterion::HistoryCriterion(const unsigned T, const unsigned ST, const OutputType HT, const double MH, uvec&& IDX)
    : Criterion(T, ST)
    , indices(IDX)
    , limit(MH)
    , history_type(HT) {}

int HistoryCriterion::process(const shared_ptr<DomainBase>& D) {
    if(indices.is_empty())
        suanpan::for_all(D->get_element_pool(), [&](const shared_ptr<Element>& t_element) {
            for(auto& I : t_element->record(history_type))
                for(const auto J : I)
                    if(check(J)) {
                        D->disable_element(t_element->get_tag());
                        return;
                    }
        });
    else
        suanpan::for_all(D->get_element_pool(), [&](const shared_ptr<Element>& t_element) {
            for(auto& I : t_element->record(history_type)) {
                uword valid_count{0};
                for(const auto J : indices)
                    if(J < I.size() && check(I(J))) valid_count += 1;
                if(valid_count == accu(indices < I.size())) {
                    D->disable_element(t_element->get_tag());
                    return;
                }
            }
        });

    return D->soft_restart();
}

bool MaxHistory::check(const double target) const { return target > limit; }

unique_ptr<Criterion> MaxHistory::unique_copy() { return std::make_unique<MaxHistory>(*this); }

bool MinHistory::check(const double target) const { return target < limit; }

unique_ptr<Criterion> MinHistory::unique_copy() { return std::make_unique<MinHistory>(*this); }
