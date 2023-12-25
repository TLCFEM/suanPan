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
 * @class LinearHardeningNM
 * @brief A LinearHardeningNM class.
 * @author tlc
 * @date 22/06/2022
 * @version 0.1.0
 * @file LinearHardeningNM.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef LINEARHARDENINGNM_H
#define LINEARHARDENINGNM_H

#include "NonlinearNM.h"

class LinearHardeningNM : public NonlinearNM {
    const double isotropic_modulus;
    const double kinematic_modulus;

    [[nodiscard]] int compute_local_integration(vec&, mat&) override;

protected:
    [[nodiscard]] vec compute_h(double) const override;
    [[nodiscard]] vec compute_dh(double) const override;

public:
    LinearHardeningNM(
        unsigned, // tag
        double,   // axial rigidity
        double,   // flexural rigidity
        double,   // isotropic hardening modulus
        double,   // kinematic hardening modulus
        double,   // linear density
        vec&&
    );
    LinearHardeningNM(
        unsigned, // tag
        double,   // axial rigidity
        double,   // flexural rigidity
        double,   // flexural rigidity
        double,   // isotropic hardening modulus
        double,   // kinematic hardening modulus
        double,   // linear density
        vec&&
    );
};

#endif

//! @}
