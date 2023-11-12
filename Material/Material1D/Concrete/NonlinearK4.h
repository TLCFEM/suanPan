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
 * @class NonlinearK4
 * @brief A ConcreteK4 material class.
 *
 * @author tlc
 * @date 05/09/2023
 * @version 0.1.0
 * @file NonlinearK4.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef NONLINEARK4_H
#define NONLINEARK4_H

#include <Material/Material1D/Material1D.h>

struct DataNonlinearK4 {
    const double elastic_modulus, hardening_k;
};

class NonlinearK4 : protected DataNonlinearK4, public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const bool apply_damage, apply_crack_closing, objective_damage;

    [[nodiscard]] virtual vec2 compute_tension_backbone(double) const = 0;
    [[nodiscard]] virtual vec2 compute_compression_backbone(double) const = 0;

    [[nodiscard]] virtual vec2 compute_tension_damage(double) const = 0;
    [[nodiscard]] virtual vec2 compute_compression_damage(double) const = 0;

    int compute_plasticity();
    void compute_crack_close_branch();

protected:
    [[nodiscard]] double objective_scale(double, double) const;

public:
    NonlinearK4(
        unsigned, // tag
        double,   // elastic modulus
        double,   // hardening parameter
        double,   // density
        bool,     // apply damage
        bool,     // apply crack closing
        bool      // objective damage
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

struct DataConcreteK4 {
    const double hardening_t, hardening_d;
    const double f_t, f_c, k_peak, f_y;
    const double zeta_t, zeta_c;
    const double hardening_c = (f_c - f_y) / k_peak;
};

class ConcreteK4 final : protected DataConcreteK4, public NonlinearK4 {
    [[nodiscard]] vec2 compute_tension_backbone(double) const override;
    [[nodiscard]] vec2 compute_compression_backbone(double) const override;

    [[nodiscard]] vec2 compute_tension_damage(double) const override;
    [[nodiscard]] vec2 compute_compression_damage(double) const override;

public:
    ConcreteK4(
        unsigned,    // tag
        double,      // elastic modulus
        double,      // hardening parameter
        vec&&,       // parameters
        double = 0., // density
        bool = true, // apply damage
        bool = true, // apply crack closing
        bool = false // objective damage
    );

    unique_ptr<Material> get_copy() override;
};

#endif

//! @}
