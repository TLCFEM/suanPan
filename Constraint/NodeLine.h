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
/**
 * @class NodeLine
 * @brief A NodeLine class.
 *
 * The `NodeLine` constraint.
 *
 * @author tlc
 * @date 02/05/2021
 * @version 0.1.0
 * @file NodeLine.h
 * @addtogroup Constraint
 * @{
 */

#ifndef NODELINE_H
#define NODELINE_H

#include "Constraint.h"

class NodeLine final : public Constraint {
    static const mat rotation;

    std::vector<vec> get_position(const shared_ptr<DomainBase>&);

public:
    NodeLine(unsigned, unsigned, unsigned, uvec&&);

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;

    void update_status(const vec&) override;
    void commit_status() override;
    void clear_status() override;
    void reset_status() override;
};

#endif

//! @}
