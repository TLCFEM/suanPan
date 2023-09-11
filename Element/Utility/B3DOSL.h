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
 * @class B3DOSL
 * @brief A B3DOSL class.
 * 
 * Order of local quantities:
 *   uniform axial
 *   strong axis bending near node
 *   strong axis bending far node
 *   weak axis bending near node
 *   weak axis bending far node
 *   torsion near node
 *   torsion far node
 *   warping near node
 *   warping far node
 *
 * @author tlc
 * @date 10/09/2023
 * @version 0.1.0
 * @file B3DOSL.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef B3DOSL_H
#define B3DOSL_H

#include "B3DL.h"

class B3DOSL final : public B3DL {
    static const span sa, sb, sc, sd;

public:
    using B3DL::B3DL;

    [[nodiscard]] unsigned input_size() const override;
    [[nodiscard]] unsigned output_size() const override;

    unique_ptr<Orientation> get_copy() override;

    [[nodiscard]] vec to_local_vec(const vec&) const override;
    [[nodiscard]] vec to_global_vec(const vec&) const override;
    [[nodiscard]] mat to_global_mass_mat(const mat&) const override;
    [[nodiscard]] mat to_global_stiffness_mat(const mat&) const override;
};

#endif

//! @}
