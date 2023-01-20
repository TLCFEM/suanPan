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
 * @class ShellBase
 * @brief A ShellBase class.
 *
 * @author tlc
 * @date 08/02/2018
 * @version 0.1.0
 * @file ShellBase.h
 * @addtogroup Shell
 * @ingroup Element
 * @{
 */

#ifndef SHELLBASE_H
#define SHELLBASE_H

#include <Element/MaterialElement.h>

class ShellBase : public MaterialElement2D {
protected:
    static const uvec m_dof, p_dof;

    static vec reshuffle(const vec&, const vec&);
    static mat reshuffle(const mat&, const mat&);

    mat trans_mat;

    void direction_cosine();
    [[nodiscard]] mat get_local_coordinate() const;
    vec& transform_from_local_to_global(vec&) const;
    vec& transform_from_global_to_local(vec&) const;
    mat& transform_from_local_to_global(mat&) const;
    vec transform_from_local_to_global(vec&&) const;
    vec transform_from_global_to_local(vec&&) const;
    mat transform_from_local_to_global(mat&&) const;

public:
    using MaterialElement2D::MaterialElement2D;
};

#endif

//! @}
