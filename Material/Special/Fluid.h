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
 * @class Fluid
 * @brief A Fluid class.
 * @author tlc
 * @date 20/11/2021
 * @file Fluid.h
 * @addtogroup Special
 * @ingroup Material
 * @{
 */

#ifndef FLUID_H
#define FLUID_H

#include <Material/Material.h>

struct DataFluid {
    const double bulk_modulus; // bulk modulus
};

class Fluid final : DataFluid, public Material {
public:
    Fluid(unsigned, // tag
          double,   // bulk modulus
          double    // density
        );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
