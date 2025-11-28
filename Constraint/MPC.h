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
/**
 * @class MPC
 * @brief A MPC class.
 *
 * @author tlc
 * @date 05/09/2017
 * @version 0.1.0
 * @file MPC.h
 * @addtogroup Constraint
 * @{
 */

#ifndef MPC_H
#define MPC_H

#include "Constraint.h"

class MPC final : public Constraint {
public:
    using Pack = std::vector<std::tuple<uword, Node::DOF, double>>;

private:
    const double magnitude;
    const Pack pool;

public:
    MPC(unsigned, unsigned, double, Pack&&);

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
