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
 * @class IsotropicNonlinearElastic3D
 * @brief The IsotropicNonlinearElastic3D class.
 *
 * @author tlc
 * @date 23/09/2020
 * @version 1.0.0
 * @file IsotropicNonlinearElastic3D.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef ISOTROPICNONLINEARELASTIC3D_H
#define ISOTROPICNONLINEARELASTIC3D_H

#include <Material/Material3D/Material3D.h>

class IsotropicNonlinearElastic3D : public Material3D {
    static constexpr double two_third = 2. / 3.;
    static const mat unit_dev_tensor;
    static const mat unit_unit;

    /**
     * \brief compute the derivatives of potential energy w.r.t. volumetric strain and equivalent strain squared
     *        (4) may be equal to (5) for most cases but not strictly always
     * \return (0) pwpm (1) pwpd (2) ppwppm (3) ppwppd (4) ppwpmpd (5) ppwpdpm
     */
    virtual vec compute_derivative(double, double) = 0;

public:
    explicit IsotropicNonlinearElastic3D(
        unsigned,   // tag
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
