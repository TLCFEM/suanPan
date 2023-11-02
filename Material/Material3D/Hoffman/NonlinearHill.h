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
 * @class NonlinearHill
 * @brief The NonlinearHill class.
 *
 * @author tlc
 * @date 20/01/2019
 * @version 0.2.0
 * @file NonlinearHill.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARHILL_H
#define NONLINEARHILL_H

#include "NonlinearHoffman.h"

class NonlinearHill : public NonlinearHoffman {
public:
    NonlinearHill(
        unsigned,   // tag
        vec&&,      // elastic modulus
        vec&&,      // poissons ratio
        vec&&,      // yield stress
        double = 0. // density
    );
};

#endif

//! @}
