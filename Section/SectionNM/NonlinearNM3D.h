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
 * @class NonlinearNM3D
 * @brief A NonlinearNM3D class.
 * @author tlc
 * @date 30/11/2021
 * @version 0.1.0
 * @file NonlinearNM3D.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef NONLINEARNM3D_H
#define NONLINEARNM3D_H

#include "SectionNM3D.h"

class NonlinearNM3D : public SectionNM3D {
    static constexpr unsigned max_iteration = 20;
    static const uvec sai, saj, sbi, sbj, sa, sb;

    const bool has_kinematic;

    const unsigned local_size = has_kinematic ? 12u : 7u;

    const uvec si{local_size - 2llu}, sj{local_size - 1llu};

    const span spa = span(0llu, local_size - 2llu);

    mat plastic_weight, kin_weight;

    [[nodiscard]] virtual double compute_h(double) const = 0;
    [[nodiscard]] virtual double compute_dh(double) const = 0;

    [[nodiscard]] virtual double compute_f(const vec&) const = 0;
    [[nodiscard]] virtual vec compute_df(const vec&) const = 0;
    [[nodiscard]] virtual mat compute_ddf(const vec&) const = 0;

protected:
    void initialize_weight(const vec&, double);

public:
    NonlinearNM3D(unsigned, // tag
                  double,
                  double,
                  double,
                  double,
                  double);

    int update_trial_status(const vec&) override;

    vector<vec> record(OutputType) override;
};

#endif

//! @}
