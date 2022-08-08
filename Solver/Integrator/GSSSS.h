/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "Integrator.h"

class GSSSS : public Integrator {
protected:
    const double L1, L2, L4;

    double L3 = 0., L5 = 0.;
    double W1 = 0., W1G1 = 0., W2G2 = 0., W3G3 = 0., W1G4 = 0., W2G5 = 0., W1G6 = 0.;

    double DT = 0.;

    double XPV2 = 0., XPV3 = 0., XPA2 = 0., XPA3 = 0., XCVD = 0., XCAD = 0.;

    // ReSharper disable once CppMemberFunctionMayBeStatic
    template<typename T> void generate_constants(double, double, double) { throw invalid_argument("need a proper scheme"); }

public:
    explicit GSSSS(unsigned);

    void assemble_resistance() override;
    void assemble_matrix() override;

    [[nodiscard]] int process_load() override;
    [[nodiscard]] int process_constraint() override;
    [[nodiscard]] int process_load_resistance() override;
    [[nodiscard]] int process_constraint_resistance() override;

    int update_trial_status() override;

    void stage_status() override;

    void update_parameter(double) override;
    void update_compatibility() const override;

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
