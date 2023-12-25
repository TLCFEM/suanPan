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
 * @class LinearSpring2D
 * @brief A LinearSpring2D class.
 *
 * @author tlc
 * @date 18/08/2022
 * @version 0.1.0
 * @file LinearSpring2D.h
 * @addtogroup Constraint
 * @{
 */

#ifndef LINEARSPRING2D_H
#define LINEARSPRING2D_H

#include "ParticleCollision2D.h"

class LinearSpring2D final : public ParticleCollision2D {
    [[nodiscard]] double compute_f(double) const override;
    [[nodiscard]] double compute_df(double) const override;

public:
    using ParticleCollision2D::ParticleCollision2D;
};

#endif

//! @}
