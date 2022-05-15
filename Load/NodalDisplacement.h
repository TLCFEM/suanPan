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
 * @class NodalDisplacement
 * @brief A NodalDisplacement class.
 *
 * The NodalDisplacement class is in charge of handling displacement load.
 *
 * @author tlc
 * @date 17/09/2017
 * @version 0.1.0
 * @file NodalDisplacement.h
 * @addtogroup Load
 * @{
 */

#ifndef NODALDISPLACEMENT_H
#define NODALDISPLACEMENT_H

#include <Load/Load.h>

class NodalDisplacement : public Load {
public:
    explicit NodalDisplacement(unsigned = 0, // tag
                               unsigned = 0, // step tag
                               double = 0.,  // magnitude
                               uvec&& = {},  // node tags
                               unsigned = 0, // dof tag
                               unsigned = 0  // amplitude tag
    );
    NodalDisplacement(unsigned,    // tag
                      unsigned,    // step tag
                      double,      // magnitude
                      uvec&&,      // node tags
                      uvec&&,      // dof tags
                      unsigned = 0 // amplitude tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;
};

#endif // DISPLACEMENTLOAD_H

//! @}
