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
 * @class T2DC
 * @brief A T2DC class.
 * @author tlc
 * @date 27/06/2018
 * @version 0.1.0
 * @file T2DC.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef T2DC_H
#define T2DC_H

#include <Element/Utility/T2DL.h>

class T2DC final : public T2DL {
protected:
    void update_transformation() override;

public:
    using T2DL::T2DL;

    [[nodiscard]] bool is_nlgeom() const override;

    unique_ptr<Orientation> get_copy() override;

    [[nodiscard]] mat to_global_geometry_mat(const mat&) const override;
};

#endif

//! @}
