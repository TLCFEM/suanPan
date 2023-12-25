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
 * @class BoucWen
 * @brief The BoucWen class.
 *
 * The BoucWen model.
 *
 * @author tlc
 * @date 24/08/2020
 * @version 1.0.0
 * @file BoucWen.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BOUCWEN_H
#define BOUCWEN_H

#include <Material/Material1D/Material1D.h>

struct DataBoucWen {
    const double elastic_modulus;
    const double yield_stress;
    const double hardening;
    const double beta;
    const double n;
};

class BoucWen final : protected DataBoucWen, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const double modulus_a = hardening * elastic_modulus;
    const double modulus_b = yield_stress - hardening * yield_stress;
    const double yield_strain = yield_stress / elastic_modulus;
    const double gamma = 1. - beta;

public:
    BoucWen(
        unsigned, // tag
        vec&&     // parameter
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
