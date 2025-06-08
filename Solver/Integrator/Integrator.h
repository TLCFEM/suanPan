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
 * @class Integrator
 * @brief The Integrator class is basically a wrapper of the DomainBase class
 * with regard to some status changing methods.
 *
 * By default, the Step object calls DomainBase(Workshop) object to update
 * displacement/resistance/stiffness independently. When it comes to dynamic
 * analysis (time integration is involved), it is necessary to compute the
 * equivalent load/stiffness by combining several quantities.
 *
 * The Integrator object is acting like an agent between Workshop and Step, that
 * can modify corresponding quantities to account for dynamic effect.
 *
 * @author tlc
 * @date 27/08/2017
 * @version 0.1.2
 * @file Integrator.h
 * @addtogroup Integrator
 * @ingroup Analysis
 * @{
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Domain/Tag.h>

class DomainBase;

enum class IntegratorType {
    Implicit,
    Explicit
};

class Integrator : public UniqueTag {
    bool time_step_switch = true;
    bool matrix_assembled_switch = false;

    std::weak_ptr<DomainBase> database;

protected:
    virtual void update_parameter(double);

public:
    explicit Integrator(unsigned = 0);

    void set_domain(const std::weak_ptr<DomainBase>&);
    [[nodiscard]] shared_ptr<DomainBase> get_domain() const;

    virtual int initialize();

    [[nodiscard]] virtual constexpr IntegratorType type() const { return IntegratorType::Implicit; }

    // ! some multistep integrators may require fixed time step for some consecutive sub-steps
    void set_time_step_switch(bool);
    [[nodiscard]] bool allow_to_change_time_step() const;

    // ! manually set switch after assembling global matrix
    void set_matrix_assembled_switch(bool);
    [[nodiscard]] bool matrix_is_assembled() const;

    [[nodiscard]] virtual bool has_corrector() const;
    [[nodiscard]] virtual bool time_independent_matrix() const;

    [[nodiscard]] virtual int process_load();
    [[nodiscard]] virtual int process_constraint();
    [[nodiscard]] int process_criterion() const;
    [[nodiscard]] int process_modifier() const;
    [[nodiscard]] virtual int process_load_resistance();
    [[nodiscard]] virtual int process_constraint_resistance();

    void record() const;

    virtual void assemble_resistance();
    virtual void assemble_matrix();

    virtual vec get_force_residual();
    virtual vec get_displacement_residual();
    [[nodiscard]] vec get_auxiliary_residual() const;
    virtual sp_mat get_reference_load();

    [[nodiscard]] virtual const vec& get_trial_displacement() const;

    void update_load() const;
    void update_constraint() const;

    void update_trial_load_factor(double) const;
    void update_trial_load_factor(const vec&) const;
    virtual void update_from_ninja();

    void update_trial_time(double);
    virtual void update_incre_time(double);

    virtual int update_trial_status(bool);
    virtual int correct_trial_status();

    int sync_status(bool);

    virtual int update_internal(const mat&);

    mat solve(const mat&);
    mat solve(const sp_mat&);
    mat solve(mat&&);
    mat solve(sp_mat&&);
    virtual int solve(mat&, const mat&);
    virtual int solve(mat&, const sp_mat&);
    virtual int solve(mat&, mat&&);
    virtual int solve(mat&, sp_mat&&);

    void erase_machine_error(vec&) const;

    void stage_and_commit_status();

    void stage_status() const;
    virtual void commit_status();
    virtual void clear_status();
    virtual void reset_status();

    virtual vec from_incre_velocity(const vec&, const uvec&);     // obtain target displacement from increment of velocity
    virtual vec from_incre_acceleration(const vec&, const uvec&); // obtain target displacement from increment of acceleration
    virtual vec from_total_velocity(const vec&, const uvec&);
    virtual vec from_total_acceleration(const vec&, const uvec&);
    vec from_incre_velocity(double, const uvec&);
    vec from_incre_acceleration(double, const uvec&);
    vec from_total_velocity(double, const uvec&);
    vec from_total_acceleration(double, const uvec&);
};

class ImplicitIntegrator : public Integrator {
public:
    using Integrator::Integrator;

    [[nodiscard]] constexpr IntegratorType type() const override { return IntegratorType::Implicit; }

    [[nodiscard]] bool time_independent_matrix() const override;
};

class ExplicitIntegrator : public Integrator {
public:
    using Integrator::Integrator;

    [[nodiscard]] constexpr IntegratorType type() const override { return IntegratorType::Explicit; }

    void assemble_resistance() override;
    void assemble_matrix() override;

    [[nodiscard]] const vec& get_trial_displacement() const override;

    void update_from_ninja() override;

    int solve(mat&, const mat&) override;
    int solve(mat&, const sp_mat&) override;
    int solve(mat&, mat&&) override;
    int solve(mat&, sp_mat&&) override;

    vec from_incre_velocity(const vec&, const uvec&) override;

    vec from_incre_acceleration(const vec&, const uvec&) override; // obtain target acceleration from increment of acceleration
    vec from_total_acceleration(const vec&, const uvec&) override;
};

#endif

//! @}
