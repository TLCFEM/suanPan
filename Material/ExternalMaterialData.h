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
 * @class ExternalMaterialData
 * @brief A ExternalMaterialData class.
 *
 * @author tlc
 * @date 15/01/2021
 * @version 0.1.0
 * @file ExternalMaterialData.h
 * @addtogroup Material
 * @{
 */

#ifndef EXTERNALMATERIALDATA_H
#define EXTERNALMATERIALDATA_H

struct ExternalMaterialData {
    unsigned size = 0;          // indicate the dimension of material
    unsigned constant_size = 0; // indicate the number of constants

    int c_strain = -1;      // current status
    int c_strain_rate = -1; // current status
    int c_stress = -1;      // current status

    int t_strain = -1;      // trial status
    int t_strain_rate = -1; // trial status
    int t_stress = -1;      // trial status

    int c_history = -1;   // current status
    int c_stiffness = -1; // stiffness matrix
    int c_damping = -1;   // damping matrix

    int t_history = -1;   // trial status
    int t_stiffness = -1; // stiffness matrix
    int t_damping = -1;   // damping matrix

    int i_history = -1;   // initial status
    int i_stiffness = -1; // stiffness matrix
    int i_damping = -1;   // damping matrix

    double density = 0.;

    double* pool = nullptr;     // stores states, should be allocated by dll
    double* constant = nullptr; // stores model constants, should be allocated by main exe
};

#endif

//! @}
