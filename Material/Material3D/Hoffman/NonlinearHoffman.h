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
 * @class NonlinearHoffman
 * @brief The NonlinearHoffman class.
 *
 * algorithm verified at 24 April 2019 by tlc
 *
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file NonlinearHoffman.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARHOFFMAN_H
#define NONLINEARHOFFMAN_H

#include <Material/Material3D/Material3D.h>

struct DataNonlinearHoffman {
    const vec modulus, ratio, yield_stress;
};

class NonlinearHoffman : protected DataNonlinearHoffman, public Material3D {
    static constexpr double two_third = 2. / 3.;
    static const double root_two_third;
    static constexpr unsigned max_iteration = 20u;
    static const uword sa;
    static const span sb;

    mat proj_a, proj_b, elastic_a;

    [[nodiscard]] virtual double compute_k(double) const = 0;
    [[nodiscard]] virtual double compute_dk(double) const = 0;

public:
    NonlinearHoffman(
        unsigned,   // tag
        vec&&,      // elastic modulus
        vec&&,      // poissons ratio
        vec&&,      // sigma
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
