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
 * @class T2DL
 * @brief A T2DL class.
 * @author tlc
 * @date 27/06/2018
 * @version 0.1.0
 * @file T2DL.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef T2DL_H
#define T2DL_H

#include <Element/Utility/Orientation.h>

class T2DL : public Orientation {
    static const span IS, JS;

protected:
    void update_transformation() override;

public:
    explicit T2DL(unsigned = 0);

    unique_ptr<Orientation> get_copy() override;

    [[nodiscard]] vec to_local_vec(const vec&) const override;
    [[nodiscard]] vec to_global_vec(const vec&) const override;
    [[nodiscard]] mat to_global_mass_mat(const mat&) const override;
    [[nodiscard]] mat to_global_stiffness_mat(const mat&) const override;
};

#endif

//! @}
