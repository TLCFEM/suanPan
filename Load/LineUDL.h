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
/**
 * @class LineUDL
 * @brief A LineUDL class.
 *
 * The LineUDL class is in charge of handling concentrated load.
 *
 * @author tlc
 * @date 19/10/2021
 * @version 0.1.0
 * @file LineUDL.h
 * @addtogroup Load
 * @{
 */

#ifndef LINEUDL_H
#define LINEUDL_H

#include <Load/Load.h>

class LineUDL : public Load {
protected:
    const uword dimension;

public:
    LineUDL(unsigned, // tag
            unsigned, // start step tag
            double,   // magnitude
            uvec&&,   // node tags
            unsigned, // dof tag
            unsigned, // amplitude tag
            uword     // dimension
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

class LineUDL2D final : public LineUDL {
public:
    LineUDL2D(unsigned, // tag
              unsigned, // start step tag
              double,   // magnitude
              uvec&&,   // node tags
              unsigned, // dof tag
              unsigned  // amplitude tag
    );

    int process(const shared_ptr<DomainBase>&) override;
};

class LineUDL3D final : public LineUDL {
public:
    LineUDL3D(unsigned, // tag
              unsigned, // start step tag
              double,   // magnitude
              uvec&&,   // node tags
              unsigned, // dof tag
              unsigned  // amplitude tag
    );

    int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
