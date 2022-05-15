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
 * @class Optimization
 * @brief A Optimization class.
 *
 * @author tlc
 * @date 12/09/2020
 * @version 0.1.0
 * @file Optimization.h
 * @addtogroup Step
 * @{
 */

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "Static.h"

class Optimization final : public Static {
public:
    using Static::Static;

    int initialize() override;

    int analyze() override;
};

#endif

//! @}
