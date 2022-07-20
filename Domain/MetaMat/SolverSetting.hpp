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

#ifndef SOLVERSETTING_HPP
#define SOLVERSETTING_HPP

#include "Preconditioner/Preconditioner.h"

enum class Precision {
    MIXED,
    FULL
};

enum class IterativeSolver {
    BICGSTAB,
    GMRES,
    NONE
};

template<sp_d data_t> struct SolverSetting {
    int restart = 20;
    int max_iteration = 200;
    data_t tolerance = std::is_same_v<data_t, float> ? 1E-7 : 1E-14;
    unsigned iterative_refinement = 5;
    Precision precision = Precision::FULL;
    IterativeSolver iterative_solver = IterativeSolver::NONE;
    unique_ptr<Preconditioner> preconditioner = nullptr;

    SolverSetting() {}

    SolverSetting(const SolverSetting& other)
        : restart(other.restart)
        , max_iteration(other.max_iteration)
        , tolerance(other.tolerance)
        , iterative_refinement(other.iterative_refinement)
        , precision(other.precision)
        , iterative_solver(other.iterative_solver)
        , preconditioner(other.preconditioner ? other.preconditioner->get_copy() : nullptr) {}

    SolverSetting(SolverSetting&&) noexcept = delete;

    SolverSetting& operator=(const SolverSetting& other) {
        if(this == &other) return *this;
        restart = other.restart;
        max_iteration = other.max_iteration;
        tolerance = other.tolerance;
        iterative_refinement = other.iterative_refinement;
        precision = other.precision;
        iterative_solver = other.iterative_solver;
        preconditioner = other.preconditioner ? other.preconditioner->get_copy() : nullptr;
        return *this;
    }

    SolverSetting& operator=(SolverSetting&&) noexcept = delete;

    ~SolverSetting() = default;
};

#endif
