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
 * @class ElementalModal
 * @brief A ElementalModal damping class.
 * @author tlc
 * @date 20/12/2018
 * @version 0.2.0
 * @file ElementalModal.h
 * @addtogroup Modifier
 * @{
 */

#ifndef ELEMENTALMODAL_H
#define ELEMENTALMODAL_H

#include <Element/Modifier/Modifier.h>

class ElementalModal final : public ModifierDynamics {
    const double cut_off_freq, damping;

public:
    ElementalModal(unsigned, double, double, uvec&& = {});

    int update_status() override;
};

#endif

//! @}
