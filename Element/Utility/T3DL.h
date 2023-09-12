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
 * @class T3DL
 * @brief A T3DL class.
 * @author tlc
 * @date 27/06/2018
 * @version 0.1.0
 * @file T3DL.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef T3DL_H
#define T3DL_H

#include <Element/Utility/Orientation.h>

class T3DL : public Orientation {
protected:
    static const span IS, JS;

    void update_transformation() override;

public:
    explicit T3DL(unsigned = 0);

    [[nodiscard]] unsigned global_size() const override;
    [[nodiscard]] unsigned local_size() const override;

    unique_ptr<Orientation> get_copy() override;

    [[nodiscard]] vec to_local_vec(const vec&) const override;
    [[nodiscard]] vec to_global_vec(const vec&) const override;
    [[nodiscard]] mat to_global_mass_mat(const mat&) const override;
    [[nodiscard]] mat to_global_stiffness_mat(const mat&) const override;
};

#endif

//! @}
