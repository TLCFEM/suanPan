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
 * @class GSSSS
 * @brief A GSSSS class defines a solver using GSSSS algorithm.
 *
 * Advances in Computational Dynamics of Particles, Materials and Structures
 *
 * ISBN: 978-1-119-96692-0
 *
 * @author tlc
 * @date 15/04/2022
 * @version 0.1.0
 * @file GSSSS.h
 * @addtogroup Integrator
 * @{
 */

#ifndef GSSSS_H
#define GSSSS_H

#include "../Integrator.h"

class GSSSS : public ImplicitIntegrator {
protected:
    void update_parameter(double) override;

    [[nodiscard]] int process_load_impl(bool) override;
    [[nodiscard]] int process_constraint_impl(bool) override;

    static constexpr double L1{1.}, L2{.5}, L4{1.};

    double L3 = 0., L5 = 0.;
    double W1 = 0., W3G3 = 0., W2G5 = 0., W1G6 = 0.;

    double DT = 0.;

    double C0{0.}, C1{0.}, C2{0.}, C3{0.}, C4{0.}, XD{0.}, XV{0.}, XA{0.};

    // ReSharper disable once CppMemberFunctionMayBeStatic
    template<typename T> void generate_constants(double, double, double) { throw std::invalid_argument("need a proper scheme"); }

public:
    explicit GSSSS(unsigned);

    void assemble_resistance() override;
    void assemble_matrix() override;

    vec get_force_residual() override;
    vec get_displacement_residual() override;
    sp_mat get_reference_load() override;

    int update_trial_status(bool) override;

    vec from_incre_velocity(const vec&, const uvec&) override;
    vec from_incre_acceleration(const vec&, const uvec&) override;

    void print() override;
};

class GSSSSU0 final : public GSSSS {
public:
    GSSSSU0(unsigned, vec&&);
};

class GSSSSV0 final : public GSSSS {
public:
    GSSSSV0(unsigned, vec&&);
};

class GSSSSOptimal final : public GSSSS {
public:
    GSSSSOptimal(unsigned, double);
};

#endif

//! @}
