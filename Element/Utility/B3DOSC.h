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
 * @class B3DOSC
 * @brief A B3DOSC class.
 * 
 * B3DOSC is a corotational transformation for 3D beam elements.
 * 
 * The implementation is mainly based on de Souza's thesis.
 * 
 * Force-based Finite Element for Large Displacement Inelastic Analysis of Frames
 * 
 * @author tlc
 * @date 16/12/2021
 * @version 0.1.0
 * @file B3DOSC.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef B3DOSC_H
#define B3DOSC_H

#include "B3DC.h"

class B3DOSC : public B3DC {
protected:
    void update_transformation() override;

    [[nodiscard]] unsigned nodal_size() const override;

public:
    B3DOSC(unsigned, vec&&);

    [[nodiscard]] OrientationType get_orientation_type() const override;

    unique_ptr<Orientation> get_copy() override;

    [[nodiscard]] vec to_local_vec(const vec&) const override;
    [[nodiscard]] mat to_global_geometry_mat(const mat&) const override;
};

#endif

//! @}
