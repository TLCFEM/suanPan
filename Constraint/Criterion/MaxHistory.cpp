/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "MaxHistory.h"
#include <Domain/DomainBase.h>
#include <Element/Element.h>

MaxHistory::MaxHistory(const unsigned T, const unsigned ST, const OutputType HT, const double MH)
    : Criterion(T, ST)
    , history_type(HT)
    , max_history(MH) {}

unique_ptr<Criterion> MaxHistory::get_copy() { return make_unique<MaxHistory>(*this); }

int MaxHistory::process(const shared_ptr<DomainBase>& D) {
    const auto& t_element_pool = D->get_element_pool();

    suanpan_for_each(t_element_pool.cbegin(), t_element_pool.cend(), [&](const shared_ptr<Element>& t_element) {
        for(auto& I : t_element->record(history_type))
            for(const auto& J : I)
                if(J > max_history) {
                    D->disable_element(t_element->get_tag());
                    return;
                }
    });

    return D->soft_restart();
}
