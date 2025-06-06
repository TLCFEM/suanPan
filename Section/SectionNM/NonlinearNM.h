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
 * @class NonlinearNM
 * @brief A NonlinearNM class.
 * @author tlc
 * @date 22/06/2022
 * @version 0.1.0
 * @file NonlinearNM.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef NONLINEARNM_H
#define NONLINEARNM_H

#include "SectionNM.h"

struct DataNonlinearNM {
    const double EA, EIS, EIW;
    const vec yield_force;
};

class NonlinearNM : protected DataNonlinearNM, public SectionNM {
    [[nodiscard]] virtual int compute_local_integration(vec&, mat&) = 0;

protected:
    static constexpr unsigned max_iteration = 20u;

    const vec yield_diag;

    const mat ti, tj;

    const bool has_kinematic;

    const unsigned n_size;                        // nodal dof size
    const unsigned d_size = 2llu * n_size - 1llu; // element dof size
    const unsigned g_size;                        // global jacobian size

    const uvec ni, nj, ga, gb, gc, gd, ge;

    [[nodiscard]] virtual vec compute_h(double) const = 0;
    [[nodiscard]] virtual vec compute_dh(double) const = 0;

    [[nodiscard]] virtual double compute_f(const vec&, const vec&) const = 0;
    [[nodiscard]] virtual vec compute_df(const vec&, const vec&) const = 0;
    [[nodiscard]] virtual mat compute_ddf(const vec&, const vec&) const = 0;

public:
    NonlinearNM(
        unsigned, // tag
        double,   // axial rigidity
        double,   // flexural rigidity
        bool,     // kinematic hardening modulus
        double,   // linear density
        vec&&
    );
    NonlinearNM(
        unsigned, // tag
        double,   // axial rigidity
        double,   // flexural rigidity
        double,   // flexural rigidity
        bool,     // kinematic hardening modulus
        double,   // linear density
        vec&&
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    std::vector<vec> record(OutputType) override;
};

#endif

//! @}
