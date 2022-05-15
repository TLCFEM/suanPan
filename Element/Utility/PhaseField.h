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
 * @class PhaseField
 * @brief A PhaseField class.
 * @author tlc
 * @date 02/12/2020
 * @version 0.1.0
 * @file PhaseField.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef PHASEFIELD_H
#define PHASEFIELD_H

#include <Material/Material.h>

struct PhaseFieldData {
    double strain_energy = 0.;
    double t_strain_energy = 0.;
    double c_strain_energy = 0.;
    double dev_strain_energy = 0.;
    double t_dev_strain_energy = 0.;
    double c_dev_strain_energy = 0.;
    double vol_strain_energy = 0.;
    double t_vol_strain_energy = 0.;
    double c_vol_strain_energy = 0.;
};

struct PhaseField : PhaseFieldData {
    void commit_status(const unique_ptr<Material>&);

    void clear_status();
};

#endif

//! @}
