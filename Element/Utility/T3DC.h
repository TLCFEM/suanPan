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
 * @class T3DC
 * @brief A T3DC class.
 * @author tlc
 * @date 27/06/2018
 * @version 0.1.0
 * @file T3DC.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef T3DC_H
#define T3DC_H

#include <Element/Utility/T3DL.h>

class T3DC final : public T3DL {
protected:
    void update_transformation() override;

public:
    explicit T3DC(unsigned = 0);

    [[nodiscard]] bool is_nlgeom() const override;

    unique_ptr<Orientation> get_copy() override;

    [[nodiscard]] mat to_global_geometry_mat(const mat&) const override;
};

#endif

//! @}
