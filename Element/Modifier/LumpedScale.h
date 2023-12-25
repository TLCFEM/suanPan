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
 * @class LumpedScale
 * @brief A LumpedScale mass class.
 * @author tlc
 * @date 20/12/2018
 * @version 0.2.0
 * @file LumpedScale.h
 * @addtogroup Modifier
 * @{
 */

#ifndef LUMPEDSCALE_H
#define LUMPEDSCALE_H

#include <Element/Modifier/Modifier.h>

class LumpedScale final : public Modifier {
public:
    using Modifier::Modifier;

    int update_status() override;
};

#endif

//! @}
