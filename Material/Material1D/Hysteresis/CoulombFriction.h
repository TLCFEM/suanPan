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
 * @class CoulombFriction
 * @brief The CoulombFriction class.
 *
 * @author tlc
 * @date 19/02/2020
 * @version 0.1.0
 * @file CoulombFriction.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef COULOMBFRICTION_H
#define COULOMBFRICTION_H

#include <Material/Material1D/Material1D.h>

struct DataCoulombFriction {
    const double friction_force;
    const double factor;
};

class CoulombFriction final : protected DataCoulombFriction, public Material1D {
public:
    CoulombFriction(
        unsigned, // tag
        double,   // maximum friction force
        double    // factor
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;
    int update_trial_status(const vec&, const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
