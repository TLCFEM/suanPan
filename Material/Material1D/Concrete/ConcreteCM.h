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
 * @class ConcreteCM
 * @brief A ConcreteCM material class.
 *
 * The ConcreteCM class represents a concrete material model based on the Chang & Mander concrete model.
 *
 * doi:10.1061/(ASCE)0733-9445(1988)114:8(1804)
 *
 * @author tlc
 * @date 29/07/2020
 * @version 0.1.0
 * @file ConcreteCM.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CONCRETECM_H
#define CONCRETECM_H

#include <Material/Material1D/Hysteresis/ComplexHysteresis.h>

class ConcreteCM final : public ComplexHysteresis {
    const double c_stress, c_strain, t_stress, t_strain;

    const double c_m, c_n, t_m, t_n;

    const bool linear_trans;

    [[nodiscard]] podarray<double> compute_compression_backbone(double) override;
    [[nodiscard]] podarray<double> compute_tension_backbone(double) override;
    [[nodiscard]] podarray<double> compute_compression_unload(double) override;
    [[nodiscard]] podarray<double> compute_tension_unload(double) override;
    [[nodiscard]] podarray<double> compute_compression_reload(double) override;
    [[nodiscard]] podarray<double> compute_tension_reload(double) override;
    [[nodiscard]] podarray<double> compute_compression_subunload(double) override;
    [[nodiscard]] podarray<double> compute_tension_subunload(double) override;
    [[nodiscard]] podarray<double> compute_transition(double, double, double, double, double, double, double) const;

    void update_compression_unload(double) override;
    void update_tension_unload(double) override;
    void update_connect();

public:
    ConcreteCM(
        unsigned,       // tag
        double,         // elastic modulus
        double,         // peak stress in negative
        double,         // crack stress in positive
        double,         // NC
        double,         // NT
        double = -2E-3, // peak strain in negative
        double = 1E-4,  // crack strain in positive
        bool = false,   // if to use linear transition
        double = 0.     // density
    );

    unique_ptr<Material> get_copy() override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    void print() override;
};

#endif

//! @}
