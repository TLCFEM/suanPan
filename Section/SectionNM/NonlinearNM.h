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
    const double EA, EIS, EIW, kinematic_modulus;
};

class NonlinearNM : protected DataNonlinearNM, public SectionNM {
    static constexpr unsigned max_iteration = 20;

    const bool has_kinematic;

    const uvec si, sj, sa, sb, sc;

    const unsigned n_size; // nodal dof size
    const unsigned j_size; // jacobian size
    const unsigned g_size = 2 * j_size;

    const vec elastic_diag, border;
    const mat ti, tj, rabbit;

    [[nodiscard]] virtual double compute_f(const vec&, double) const = 0;
    [[nodiscard]] virtual double compute_dh(const vec&, double) const =0;
    [[nodiscard]] virtual vec compute_df(const vec&, double) const = 0;
    [[nodiscard]] virtual mat compute_ddf(const vec&, double) const = 0;

    bool update_nodal_quantity(mat&, vec&, double, const vec&, const vec&, double, const vec&, const vec&) const;

public:
    NonlinearNM(unsigned, // tag
                double,   // axial rigidity
                double,   // flexural rigidity
                double,   // kinematic hardening modulus
                double    // linear density
    );
    NonlinearNM(unsigned, // tag
                double,   // axial rigidity
                double,   // flexural rigidity
                double,   // flexural rigidity
                double,   // kinematic hardening modulus
                double    // linear density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    vector<vec> record(OutputType) override;
};

#endif

//! @}
