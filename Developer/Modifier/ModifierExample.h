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
 * @class ModifierExample
 * @brief A ModifierExample damping class.
 * @author tlc
 * @date 20/2/2019
 * @version 0.1.0
 * @file ModifierExample.h
 * @addtogroup Modifier
 * @{
 */

#ifndef MODIFIEREXAMPLE_H
#define MODIFIEREXAMPLE_H

#include <Element/Modifier/Modifier.h>

class ModifierExample final : public Modifier {
    const double a, b;

public:
    ModifierExample(unsigned, double, double, uvec&& = {});

    int update_status() override;
};

SUANPAN_EXPORT void new_modifierexample(unique_ptr<Modifier>&, std::istringstream&);

#endif

//! @}
