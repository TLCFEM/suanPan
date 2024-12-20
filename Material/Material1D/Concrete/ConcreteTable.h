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
 * @class ConcreteTable
 * @brief A ConcreteTable material class.
 *
 * @author tlc
 * @date 08/09/2019
 * @version 0.1.0
 * @file ConcreteTable.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CONCRETETABLE_H
#define CONCRETETABLE_H

#include <Material/Material1D/Hysteresis/SimpleHysteresis.h>

class ConcreteTable final : public SimpleHysteresis {
    mat c_table, t_table;

    const double c_strain, t_strain;

    [[nodiscard]] podarray<double> compute_compression_initial_reverse() const override;
    [[nodiscard]] podarray<double> compute_tension_initial_reverse() const override;
    [[nodiscard]] podarray<double> compute_compression_backbone(double) const override;
    [[nodiscard]] podarray<double> compute_tension_backbone(double) const override;
    [[nodiscard]] double compute_compression_residual(double, double) const override;
    [[nodiscard]] double compute_tension_residual(double, double) const override;

public:
    ConcreteTable(
        unsigned,   // tag
        mat&&,      // compression table
        mat&&,      // tension table
        double,     // middle point
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
