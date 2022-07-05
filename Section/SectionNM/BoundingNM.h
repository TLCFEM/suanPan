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
 * @class BoundingNM
 * @brief A BoundingNM class.
 * @author tlc
 * @date 04/07/2022
 * @version 0.1.0
 * @file BoundingNM.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef BOUNDINGNM_H
#define BOUNDINGNM_H

#include "NonlinearNM.h"

class BoundingNM : public NonlinearNM {
    const double iso_modulus, iso_saturation, iso_decay;
    const double bounding_saturation, bounding_rate;

    [[nodiscard]] int compute_local_integration(vec&, mat&, bool, bool) override;

protected:
    [[nodiscard]] vec compute_h(double) const override;
    [[nodiscard]] vec compute_dh(double) const override;

public:
    BoundingNM(unsigned, // tag
               double,   // axial rigidity
               double,   // flexural rigidity
               double,   // isotropic modulus
               double,   // isotropic saturation
               double,   // isotropic decay
               double,   // bounding saturation
               double,   // bounding rate
               double,   // linear density
               vec&&
    );
    BoundingNM(unsigned, // tag
               double,   // axial rigidity
               double,   // flexural rigidity
               double,   // flexural rigidity
               double,   // isotropic modulus
               double,   // isotropic saturation
               double,   // isotropic decay
               double,   // bounding saturation
               double,   // bounding rate
               double,   // linear density
               vec&&
    );

};

#endif

//! @}
