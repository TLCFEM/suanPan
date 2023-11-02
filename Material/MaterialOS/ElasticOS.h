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
 * @class ElasticOS
 * @brief The ElasticOS class defines a isotropic elastic material for open section
 * problems.
 *
 * The Young's modulus is stored in `elastic_modulus`. The Poisson's ratio is
 * stored in `poissons_ratio`.
 *
 *
 * @author tlc
 * @date 15/09/2023
 * @version 1.0.0
 * @file ElasticOS.h
 * @addtogroup Material-OS
 * @{
 */

#ifndef ELASTICOS_H
#define ELASTICOS_H

#include <Material/MaterialOS/MaterialOS.h>

struct DataElasticOS {
    double elastic_modulus; // elastic modulus
    double poissons_ratio;  // poissons ratio
};

class ElasticOS final : public DataElasticOS, public MaterialOS {
public:
    ElasticOS(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // poissons ratio
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
