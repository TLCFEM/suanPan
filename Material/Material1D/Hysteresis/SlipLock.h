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
 * @class SlipLock
 * @brief The SlipLock class.
 *
 * The SlipLock model based on Menegotto--Pinto equation.
 *
 * validated @ 29/08/2019 by tlc
 *
 * @author tlc
 * @date 29/08/2019
 * @version 1.0.0
 * @file SlipLock.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef SLIPLOCK_H
#define SLIPLOCK_H

#include <Material/Material1D/Material1D.h>

struct DataSlipLock {
    const double elastic_modulus; // elastic modulus
    const double hardening_ratio; // hardening ratio
    const double yield_strain;    // yield stress
    const double R0;              // model parameters
};

class SlipLock final : DataSlipLock, public Material1D {
    static constexpr unsigned max_iteration = 20;

    const double yield_stress = yield_strain * elastic_modulus; // yield strain
public:
    SlipLock(unsigned,     // tag
             double,       // elastic modulus
             double,       // yield strain
             double = .05, // hardening ratio
             double = 20., // R0
             double = 0.   // density
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
