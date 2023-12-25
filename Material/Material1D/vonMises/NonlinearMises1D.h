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
 * @class NonlinearMises1D
 * @brief A NonlinearMises1D material class.
 * @author tlc
 * @date 21/01/2019
 * @version 0.1.0
 * @file NonlinearMises1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef NONLINEARMISES1D_H
#define NONLINEARMISES1D_H

#include <Material/Material1D/Material1D.h>

struct DataMises1D {
    const double elastic_modulus; // elastic modulus
};

class NonlinearMises1D : protected DataMises1D, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    [[nodiscard]] virtual double compute_k(double) const = 0;
    [[nodiscard]] virtual double compute_dk(double) const = 0;
    [[nodiscard]] virtual double compute_h(double) const = 0;
    [[nodiscard]] virtual double compute_dh(double) const = 0;

public:
    NonlinearMises1D(
        unsigned,   // tag
        double,     // elastic modulus
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
