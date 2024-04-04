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
 * @class DuncanSelig
 * @brief A DuncanSelig material class.
 *
 * For plane strain soil.
 *
 * @author tlc
 * @date 01/04/24
 * @version 1.0.0
 * @file DuncanSelig.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef DUNCANSELIG_H
#define DUNCANSELIG_H

#include <Material/Material2D/Material2D.h>

class DuncanSelig final : public Material2D {
    static constexpr unsigned max_iteration = 20u;

    static double dev(const vec&);
    static rowvec3 der_dev(const vec&);
    static mat33 compute_stiffness(double, double);

    double p_atm = 14.7;
    double ref_elastic = 400. * p_atm, n = .6;
    double ref_bulk = 300. * p_atm, m = .2;
    double ini_phi = .7, ten_fold_phi_diff = .1, r_f = .7, cohesion = .5;

    int project_to_surface(double&);

    std::tuple<double, double, double, rowvec3, rowvec3> compute_moduli();

public:
    DuncanSelig(
        unsigned,   // tag
        const vec&, // parameters
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
