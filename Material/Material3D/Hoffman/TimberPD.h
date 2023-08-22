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
 * @class TimberPD
 * @brief The TimberPD class.
 *
 * @author tlc
 * @date 24/08/2023
 * @version 1.0.0
 * @file TimberPD.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef TIMBERPD_H
#define TIMBERPD_H

#include <Material/Material3D/Material3D.h>

struct DataTimberPD {
    const vec modulus, ratio, yield_stress;

    mat proj_a;
    vec proj_b;

    mat hill_t, hill_c;

    double h, ini_r_t, n_t, b_t, ini_r_c, m_c, b_c;
};

class TimberPD final : protected DataTimberPD, public Material3D {
    static const unsigned max_iteration;

    static const uword sa;
    static const span sb;

    [[nodiscard]] double compute_damage_c(double) const;
    [[nodiscard]] double compute_damage_t(double) const;
    [[nodiscard]] double update_damage_t(const vec&, mat&);
    [[nodiscard]] double update_damage_c(const vec&, mat&);

public:
    TimberPD(unsigned,   // tag
             vec&&,      // elastic modulus
             vec&&,      // poissons ratio
             vec&&,      // sigma
             vec&&,      // hardening
             double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
