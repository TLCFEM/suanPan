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
 * @class Substepping
 * @brief A Substepping material class.
 * @author tlc
 * @date 12/09/2020
 * @version 0.1.0
 * @file Substepping.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef SUBSTEPPING_H
#define SUBSTEPPING_H

#include <Material/Material.h>

class Substepping final : public Material {
    const unsigned max_iteration;

    const unsigned mat_tag;

    unique_ptr<Material> trial_mat_obj, current_mat_obj;

public:
    Substepping(unsigned, // tag
                unsigned, // mat tag
                unsigned  // max iteration
    );
    Substepping(const Substepping&);
    Substepping(Substepping&&) = delete;
    Substepping& operator=(const Substepping&) = delete;
    Substepping& operator=(Substepping&&) = delete;
    ~Substepping() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    [[nodiscard]] const mat& get_initial_damping() const override;
    [[nodiscard]] const vec& get_initial_history() const override;
    [[nodiscard]] const mat& get_initial_stiffness() const override;

    const mat& get_current_damping() override;
    const mat& get_current_secant() override;
    const mat& get_current_stiffness() override;

    const mat& get_trial_damping() override;
    const mat& get_trial_secant() override;
    const mat& get_trial_stiffness() override;

    const vec& get_current_strain() override;
    const vec& get_current_strain_rate() override;
    const vec& get_current_stress() override;

    const vec& get_trial_strain() override;
    const vec& get_trial_strain_rate() override;
    const vec& get_trial_stress() override;

    void set_initial_history(const vec&) override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
