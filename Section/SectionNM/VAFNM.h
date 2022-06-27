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
 * @class VAFNM
 * @brief A VAFNM class.
 * @author tlc
 * @date 22/06/2022
 * @version 0.1.0
 * @file VAFNM.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef VAFNM_H
#define VAFNM_H

#include "NonlinearNM.h"

class VAFNM : public NonlinearNM {
    const double iso_modulus;
    const double kin_modulus;

    const double iso_saturation, iso_decay;
    const double kin_base;

    bool update_nodal_quantity(mat&, vec&, double, const vec&, const vec&, double, const vec&, const vec&) const override;

protected:
    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

public:
    VAFNM(unsigned, // tag
          double,   // axial rigidity
          double,   // flexural rigidity
          double,   // isotropic modulus
          double,   // isotropic saturation
          double,   // isotropic decay
          double,   // kinematic modulus
          double,   // kinematic base
          double,   // linear density
          vec&&
    );
    VAFNM(unsigned, // tag
          double,   // axial rigidity
          double,   // flexural rigidity
          double,   // flexural rigidity
          double,   // isotropic modulus
          double,   // isotropic saturation
          double,   // isotropic decay
          double,   // kinematic modulus
          double,   // kinematic base
          double,   // linear density
          vec&&
    );

};

#endif

//! @}
