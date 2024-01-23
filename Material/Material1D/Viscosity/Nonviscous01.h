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
 * @class Nonviscous01
 * @brief A 1D Viscosity class.
 *
 * Reference: 10.1016/j.ymssp.2024.111156
 *
 * @author tlc
 * @date 28/03/2023
 * @version 0.2.0
 * @file Nonviscous01.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef NONVISCOUS01_H
#define NONVISCOUS01_H

#include <Material/Material1D/Material1D.h>

struct DataNonviscous01 {
    const cx_vec m, s;
};

class Nonviscous01 final : protected DataNonviscous01, public Material1D {
    const double* incre_time = nullptr;

    cx_vec complex_damping;

    double accu_para{0.};
    cx_vec s_para, m_para;

public:
    Nonviscous01(
        unsigned, // tag
        cx_vec&&, // m
        cx_vec&&  // s
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;
    int update_trial_status(const vec&, const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
