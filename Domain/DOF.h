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
 * @author tlc
 * @date 20/05/2022
 * @version 0.1.0
 * @file DOF.h
 * @addtogroup Node
 * @{
 */

#ifndef DOF_H
#define DOF_H

enum class DOF : unsigned short {
    NONE,
    X,  // displacement in x direction
    Y,  // displacement in y direction
    Z,  // displacement in z direction
    RX, // rotation in x direction
    RY, // rotation in y direction
    RZ, // rotation in z direction
    DMG, //damage
    P,  // pressure
    T   // temperature
};

#endif // DOF_H

//! @}
