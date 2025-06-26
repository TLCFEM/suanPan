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
 * @class ParticleCollision3D
 * @brief A ParticleCollision3D class.
 *
 * @author tlc
 * @date 15/07/2020
 * @version 0.1.0
 * @file ParticleCollision3D.h
 * @addtogroup Constraint
 * @{
 */

#ifndef PARTICLECOLLISION3D_H
#define PARTICLECOLLISION3D_H

#include "ParticleCollision.h"

class ParticleCollision3D final : public ParticleCollision {
    struct CellList {
        int x = 0, y = 0, z = 0;
        unsigned tag = 0;
    };

    std::vector<CellList> list;

    const double space = 1.;
    const double alpha = 1.;

    [[nodiscard]] double compute_f(double) const override;
    [[nodiscard]] double compute_df(double) const override;

    int process_meta(const shared_ptr<DomainBase>&, bool) override;

public:
    ParticleCollision3D(unsigned, double, double);
};

#endif

//! @}
