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
 * @class ConcreteExp
 * @brief A ConcreteExp material class.
 *
 * @author tlc
 * @date 28/10/2019
 * @version 0.1.0
 * @file ConcreteExp.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CONCRETEEXP_H
#define CONCRETEEXP_H

#include <Material/Material1D/Hysteresis/SimpleHysteresis.h>

struct DataConcreteExp {
    const double elastic_modulus;
    const double f_t, f_c, a_t, a_c, b_t, b_c;
};

class ConcreteExp final : protected DataConcreteExp, public SimpleHysteresis {
    static constexpr unsigned max_iteration = 20u;

    [[nodiscard]] podarray<double> compute_compression_initial_reverse() const override;
    [[nodiscard]] podarray<double> compute_tension_initial_reverse() const override;
    [[nodiscard]] podarray<double> compute_compression_backbone(double) const override;
    [[nodiscard]] podarray<double> compute_tension_backbone(double) const override;
    [[nodiscard]] double compute_compression_residual(double, double) const override;
    [[nodiscard]] double compute_tension_residual(double, double) const override;

public:
    ConcreteExp(unsigned,   // tag
                double,     // elastic modulus
                double,     // f_t
                double,     // a_t
                double,     // g_t
                double,     // f_c
                double,     // a_c
                double,     // g_c
                double,     // middle point
                double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
