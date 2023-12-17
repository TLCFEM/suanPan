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
/**
 * @class Embed
 * @brief A Embed class.
 *
 * @author tlc
 * @date 12/08/2020
 * @version 0.1.0
 * @file Embed.h
 * @addtogroup Constraint
 * @{
 */

#ifndef EMBED2D_H
#define EMBED2D_H

#include "Constraint.h"

class Embed : public Constraint {
    static constexpr unsigned max_iteration = 20u;

    const unsigned element_tag;

public:
    Embed(
        unsigned, // unique constraint tag
        unsigned, // step tag
        unsigned, // element tag
        unsigned, // node tag
        unsigned  // dimension
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;
};

class Embed2D final : public Embed {
public:
    Embed2D(unsigned, unsigned, unsigned, unsigned);
};

class Embed3D final : public Embed {
public:
    Embed3D(unsigned, unsigned, unsigned, unsigned);
};

#endif

//! @}
