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
 * @class GroupBodyForce
 * @brief A GroupBodyForce class.
 *
 * The GroupBodyForce class is in charge of handling body force.
 *
 * @author tlc
 * @date 11/07/2020
 * @version 0.1.0
 * @file GroupBodyForce.h
 * @addtogroup Load
 * @{
 */

#ifndef GROUPBODYFORCE_H
#define GROUPBODYFORCE_H

#include "BodyForce.h"

class GroupBodyForce final : public BodyForce {
    const uvec groups;

    void update_element_tag(const shared_ptr<DomainBase>&);

public:
    explicit GroupBodyForce(unsigned = 0, // tag
                            unsigned = 0, // start step tag
                            double = 0.,  // magnitude
                            uvec&& = {},  // element tags
                            unsigned = 0, // dof tag
                            unsigned = 0  // amplitude tag
    );
    GroupBodyForce(unsigned,    // tag
                   unsigned,    // start step tag
                   double,      // magnitude
                   uvec&&,      // element tags
                   uvec&&,      // dof tags
                   unsigned = 0 // amplitude tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
