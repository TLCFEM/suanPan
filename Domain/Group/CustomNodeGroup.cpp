/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "CustomNodeGroup.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>

CustomNodeGroup::CustomNodeGroup(const unsigned T, const unsigned ET)
    : Group(T)
    , expression_tag(ET) {}

void CustomNodeGroup::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(expression_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", expression_tag);
        return;
    }

    expression = D->get_expression(expression_tag);

    if(expression->output_size() != 1) {
        suanpan_error("The assigned expression {} does not generate single output.\n", expression_tag);
        return;
    }

    std::vector<uword> pond;

    const auto input_size = expression->input_size();

    for(auto& I : D->get_node_pool()) if(as_scalar(expression->evaluate(resize(I->get_coordinate(), input_size, 1llu))) > 0.5) pond.emplace_back(I->get_tag());

    suanpan_sort(pond.begin(), pond.end());

    pool = pond;
}
