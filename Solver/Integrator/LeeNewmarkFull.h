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
 * @class LeeNewmarkFull
 * @brief A LeeNewmarkFull class defines a solver using Newmark algorithm with Lee damping model.
 *
 * doi: 10.1016/j.compstruc.2020.106423
 * doi: 10.1016/j.compstruc.2021.106663
 *
 * @author tlc
 * @date 26/11/2021
 * @version 0.1.0
 * @file LeeNewmarkFull.h
 * @addtogroup Integrator
 * @{
 */

#ifndef LEENEWMARKFULL_H
#define LEENEWMARKFULL_H

#include "LeeNewmarkBase.h"

class LeeNewmarkFull final : public LeeNewmarkBase {
public:
    enum class Type {
        T0,
        T1,
        T2,
        T3,
        T4
    };

    struct Mode {
        Type t;
        vec p;
        double zeta, omega;
    };

private:
    std::vector<Mode> damping_mode;

    // a temporary matrix holding first block nonzero entries from damping matrix
    const triplet_form<double, uword> rabbit;
    const triplet_form<double, uword> current_stiffness;
    const triplet_form<double, uword> current_mass;

    using index_tm = decltype(current_mass)::index_type;
    using index_ts = decltype(current_stiffness)::index_type;

    [[nodiscard]] uword get_amplifier() const;
    [[nodiscard]] uword get_total_size() const override;

    void update_stiffness() const override;
    void update_residual() const override;

    void formulate_block(uword&, double, double, int) const;
    void formulate_block(uword&, const std::vector<double>&, const std::vector<double>&, const std::vector<int>&) const;
    void assemble_by_mode_zero(uword&, double, double) const;
    void assemble_by_mode_one(uword&, double, double, int) const;
    void assemble_by_mode_two(uword&, double, double, int, int) const;
    void assemble_by_mode_three(uword&, double, double, double) const;
    void assemble_by_mode_four(uword&, double, double, int, int, int, int, double) const;

public:
    explicit LeeNewmarkFull(unsigned, std::vector<Mode>&&, double, double, StiffnessType);

    int initialize() override;

    int process_constraint() override;

    void print() override;
};

#endif

//! @}
