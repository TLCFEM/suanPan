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
 * @class Solver
 * @brief A Solver class defines solvers used in analysis.
 *
 * @author tlc
 * @date 27/07/2017
 * @version 0.2.1
 * @file Solver.h
 * @addtogroup Solver
 * @ingroup Analysis
 * @{
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <Domain/Tag.h>

class Converger;
class Integrator;

class Solver : public UniqueTag {
    shared_ptr<Converger> converger = nullptr;
    shared_ptr<Integrator> modifier = nullptr;

    double step_amplifier = 1.0;

protected:
    [[nodiscard]] bool constant_matrix() const;

public:
    explicit Solver(unsigned = 0);

    virtual int initialize();

    virtual int analyze() = 0;

    virtual void set_step_size(double) {}

    void set_step_amplifier(double);
    [[nodiscard]] double get_step_amplifier() const;

    void set_converger(const shared_ptr<Converger>&);
    [[nodiscard]] const shared_ptr<Converger>& get_converger() const;

    void set_integrator(const shared_ptr<Integrator>&);
    [[nodiscard]] const shared_ptr<Integrator>& get_integrator() const;
};

#endif

//! @}
