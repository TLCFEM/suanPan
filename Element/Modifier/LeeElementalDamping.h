/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class LeeElementalDamping
 * @brief A LeeElementalDamping damping class.
 * @author tlc
 * @date 20/10/2022
 * @version 0.1.0
 * @file LeeElementalDamping.h
 * @addtogroup Modifier
 * @{
 */

#ifndef LEEELEMENTALDAMPING_H
#define LEEELEMENTALDAMPING_H

#include <Element/Modifier/Modifier.h>

class LeeElementalDamping final : public Modifier {
    const double default_damping_ratio = .02;
    const double damping_ratio;

public:
    LeeElementalDamping(unsigned, double, uvec&& = {});

    int update_status() override;
};

#endif

//! @}
