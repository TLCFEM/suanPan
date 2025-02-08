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
 * @class BWBN
 * @brief The BWBN class.
 *
 * The BWBN model.
 *
 * @author tlc
 * @date 24/08/2020
 * @version 1.0.0
 * @file BWBN.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BWBN_H
#define BWBN_H

#include <Material/Material1D/Material1D.h>

struct DataBWBN {
    const vec pool;
};

class BWBN final : protected DataBWBN, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const double& elastic_modulus = pool(0);
    const double& yield_stress = pool(1);
    const double& hardening = pool(2);
    const double& beta = pool(3);
    const double& n = pool(4);
    const double& nu_i = pool(5);
    const double& nu_rate = pool(6);
    const double& eta_i = pool(7);
    const double& eta_rate = pool(8);
    const double& phi_i = pool(9);
    const double& phi_rate = pool(10);
    const double& zeta = pool(11);
    const double& a_rate = pool(12);
    const double& p = pool(13);
    const double& q = pool(14);
    const double& lambda = pool(15);

    const double modulus_a = hardening * elastic_modulus;
    const double modulus_b = yield_stress - hardening * yield_stress;
    const double modulus_c = modulus_b / elastic_modulus;
    const double yield_strain = yield_stress / elastic_modulus;

public:
    BWBN(
        unsigned, // tag
        vec&&,    // parameter
        double    // density
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
