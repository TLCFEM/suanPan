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
 * @class LeeNewmark
 * @brief A LeeNewmark class defines a solver using Newmark algorithm with Lee damping model.
 *
 * Remarks:
 *   1. Both dense and sparse storage schemes for original matrices are considered.
 *      The matrix--vector product may potentially be faster with dense formulation.
 *
 * @author tlc
 * @date 25/05/2020
 * @version 0.1.1
 * @file LeeNewmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef LEENEWMARK_H
#define LEENEWMARK_H

#include "LeeNewmarkBase.h"

class LeeNewmark : public LeeNewmarkBase {
    const vec mass_coef, stiffness_coef;

    const uword n_damping = mass_coef.n_elem;

    const double CM;

    [[nodiscard]] uword get_total_size() const override;

    void update_stiffness() const override;
    void update_residual() const override;

    virtual void initialize_mass(const shared_ptr<DomainBase>&);
    virtual void initialize_stiffness(const shared_ptr<DomainBase>&);

protected:
    shared_ptr<MetaMat<double>> current_mass = nullptr;
    shared_ptr<MetaMat<double>> current_stiffness = nullptr;
    shared_ptr<MetaMat<double>> current_geometry = nullptr;

public:
    explicit LeeNewmark(unsigned, vec&&, vec&&, double, double);

    int initialize() override;

    int process_constraint() override;
    int process_constraint_resistance() override;

    void assemble_resistance() override;

    void print() override;
};

class LeeElementalNewmark final : public LeeNewmark {
    void initialize_mass(const shared_ptr<DomainBase>&) override;
    void initialize_stiffness(const shared_ptr<DomainBase>&) override;

public:
    using LeeNewmark::LeeNewmark;
};

#endif

//! @}
