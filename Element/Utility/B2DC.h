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
 * @class B2DC
 * @brief A B2DC class.
 * @author tlc
 * @date 19/12/2021
 * @version 0.1.0
 * @file B2DC.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef B2DC_H
#define B2DC_H

#include "B2DL.h"

class B2DC final : public B2DL {
    vec original_position;

    void update_transformation() override;

public:
    using B2DL::B2DL;

    [[nodiscard]] bool is_nlgeom() const override;

    unique_ptr<Orientation> get_copy() override;

    [[nodiscard]] vec to_local_vec(const vec&) const override;
    [[nodiscard]] mat to_global_geometry_mat(const mat&) const override;
};

#endif

//! @}
