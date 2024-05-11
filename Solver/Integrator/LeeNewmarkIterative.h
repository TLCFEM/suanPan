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
 * @class LeeNewmarkIterative
 * @brief A LeeNewmarkIterative class defines a solver using Newmark algorithm with Lee damping model.
 *
 * doi:10.1016/j.compstruc.2020.106423
 * doi:10.1016/j.compstruc.2021.106663
 *
 * Remarks:
 *   1. An iterative algorithm is used. The quadratic convergence is lost.
 *
 * @author tlc
 * @date 15/01/2024
 * @version 0.1.0
 * @file LeeNewmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef LEENEWMARKITERATIVE_H
#define LEENEWMARKITERATIVE_H

#include "Newmark.h"
#include <Domain/MetaMat/MetaMat.hpp>
#include <Domain/Factory.hpp>

class LeeNewmarkIterative final : public Newmark {
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
    unsigned n_block{0};

    std::vector<Mode> damping_mode;

    shared_ptr<Factory<double>> factory = nullptr;

    shared_ptr<MetaMat<double>> current_mass = nullptr;
    shared_ptr<MetaMat<double>> current_stiffness = nullptr;

    unique_ptr<MetaMat<double>> worker = nullptr;

    void init_worker(unsigned, unsigned);

    void assemble(const shared_ptr<MetaMat<double>>&, uword, uword, double) const;

    void assemble_mass(uword, uword, double) const;
    void assemble_stiffness(uword, uword, double) const;
    void assemble_mass(const std::vector<sword>&, const std::vector<sword>&, const std::vector<double>&) const;
    void assemble_stiffness(const std::vector<sword>&, const std::vector<sword>&, const std::vector<double>&) const;

    void formulate_block(sword&, double, double, int) const;
    void formulate_block(sword&, const std::vector<double>&, const std::vector<double>&, const std::vector<int>&) const;

    [[nodiscard]] vec update_by_mode_one(double, double, int) const;
    [[nodiscard]] vec update_by_mode_two(double, double, int, int);
    [[nodiscard]] vec update_by_mode_three(double, double, double);
    [[nodiscard]] vec update_by_mode_four(double, double, int, int, int, int, double);

    void update_damping_force();

public:
    LeeNewmarkIterative(unsigned, std::vector<Mode>&&, double, double);

    int initialize() override;

    [[nodiscard]] int process_constraint() override;
    [[nodiscard]] int process_constraint_resistance() override;

    void assemble_matrix() override;

    void print() override;
};

#endif

//! @}
