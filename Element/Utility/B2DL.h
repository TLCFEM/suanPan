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
 * @class B2DL
 * @brief A B2DL class.
 * @author tlc
 * @date 27/06/2018
 * @version 0.1.0
 * @file B2DL.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef B2DL_H
#define B2DL_H

#include <Element/Utility/Orientation.h>

class B2DL : public Orientation {
protected:
    mat trans_mat;

    void form_trans_mat(const vec&);
    void update_transformation() override;

public:
    using Orientation::Orientation;

    [[nodiscard]] OrientationType get_orientation_type() const override;

    unique_ptr<Orientation> get_copy() override;

    [[nodiscard]] vec to_local_vec(const vec&) const override;
    [[nodiscard]] vec to_global_vec(const vec&) const override;
    [[nodiscard]] mat to_global_mass_mat(const mat&) const override;
    [[nodiscard]] mat to_global_stiffness_mat(const mat&) const override;
};

#endif

//! @}
