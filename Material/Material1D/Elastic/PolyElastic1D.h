/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
 * @class PolyElastic1D
 * @brief A PolyElastic1D material class.
 * @author tlc
 * @date 25/08/2020
 * @version 0.1.0
 * @file PolyElastic1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef POLYELASTIC1D_H
#define POLYELASTIC1D_H

#include <Material/Material1D/Material1D.h>

struct DataPolyElastic1D {
    const vec pool;
};

class PolyElastic1D final : protected DataPolyElastic1D, public Material1D {
public:
    PolyElastic1D(
        unsigned,   // tag
        vec&&,      // parameters
        double = 0. // density
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
