/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "ElementGroup.h"
#include <Domain/DomainBase.h>
#include <Element/Element.h>

ElementGroup::ElementGroup(const unsigned T, uvec&& R)
    : Group(T, std::forward<uvec>(R)) {}

void ElementGroup::initialize(const shared_ptr<DomainBase>& D) {
    if(!pool.empty()) return;

    const auto& e_pool = D->get_element_pool();
    pool.set_size(e_pool.size());
    for(uword I = 0; I < pool.n_elem; ++I) pool(I) = e_pool[I]->get_tag();
}

void ElementGroup::print() {
    suanpan_info("An element group.\n");
}
