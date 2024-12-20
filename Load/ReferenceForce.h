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
 * @class ReferenceForce
 * @brief A ReferenceForce class.
 *
 * The ReferenceForce class is in charge of handling concentrated reference load.
 *
 * @author tlc
 * @date 23/10/2023
 * @version 0.1.0
 * @file ReferenceForce.h
 * @addtogroup Load
 * @{
 */

#ifndef REFERENCEFORCE_H
#define REFERENCEFORCE_H

#include <Load/Load.h>

class ReferenceForce final : public Load {
public:
    ReferenceForce(
        unsigned, // tag
        unsigned, // start step tag
        double,   // magnitude
        uvec&&,   // node tags
        unsigned  // dof tag
    );

    int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
