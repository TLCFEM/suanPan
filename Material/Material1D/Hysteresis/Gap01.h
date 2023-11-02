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
 * @class Gap01
 * @brief The Gap01 class defines a gap material. This material does not
 * response to compression, but show bilinear elastic-perfectly plastic response
 * with the elastic modulus of \f$E\f$ and the yield stress of \f$\sigma_y\f$.
 * The material also shows a gap feature during reloading stage, i.e., from zero
 * stress to tension stress.
 *
 * @author tlc
 * @date 05/07/2018
 * @version 0.2.0
 * @file Gap01.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef GAP01_H
#define GAP01_H

#include <Material/Material1D/Material1D.h>

struct DataGap01 {
    const double elastic_modulus;
    const double yield_stress;
    const double gap_strain;
};

class Gap01 final : protected DataGap01, public Material1D {
public:
    Gap01(
        unsigned,    // tag
        double,      // elastic modulus
        double,      // yield stress
        double = 0., // gap strain
        double = 0.  // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
