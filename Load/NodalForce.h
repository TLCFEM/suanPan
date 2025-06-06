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
 * @class NodalForce
 * @brief A NodalForce class.
 *
 * The NodalForce class is in charge of handling concentrated load.
 *
 * @author tlc
 * @date 23/07/2017
 * @version 0.1.0
 * @file NodalForce.h
 * @addtogroup Load
 * @{
 */

#ifndef NODALFORCE_H
#define NODALFORCE_H

#include <Load/Load.h>

class NodalForce : public Load {
public:
    NodalForce(
        unsigned, // tag
        double,   // magnitude
        uvec&&,   // node tags
        uvec&&,   // dof tags
        unsigned  // amplitude tag
    );

    int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
