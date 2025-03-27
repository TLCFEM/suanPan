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
 * @class Step
 * @brief A Step class.
 * @author tlc
 * @date 27/08/2017
 * @version 0.2.1
 * @file Step.h
 * @addtogroup Step
 * @{
 */

#ifndef STEP_H
#define STEP_H

#include <Domain/Factory.hpp>
#include <Domain/Tag.h>

class DomainBase;
class Solver;
class Converger;
class Integrator;

class Step : public Tag {
    double time_period = 1.0; // time period

    double time_left = time_period;

    double max_step_size = time_period;                         // maximum step size
    double min_step_size = time_period > 0. ? 1E-8 : 0.;        // minimum step size
    double ini_step_size = time_period > 0. ? time_period : 1.; // initial step size

    unsigned max_substep = 1000; // maximum increment number

    bool fixed_step_size = false; // auto-stepping

    unsigned solver_tag = 0;
    unsigned converger_tag = 0;
    unsigned integrator_tag = 0;

protected:
    const bool symm_mat = false;
    const bool band_mat = true;
    const bool sparse_mat = false;

    SolverType system_solver = SolverType::LAPACK;
    SolverSetting<double> system_setting{};

    SolverType sub_system_solver = SolverType::LAPACK;

#ifdef SUANPAN_MAGMA
    magma_dopts magma_setting{};
#endif

    weak_ptr<DomainBase> database;
    shared_ptr<Factory<double>> factory;
    shared_ptr<Solver> solver;
    shared_ptr<Converger> tester;
    shared_ptr<Integrator> modifier;

    void configure_storage_scheme() const;

public:
    explicit Step(unsigned = 0, double = 1.);
    Step(const Step&) = delete;
    Step(Step&&) noexcept = delete;
    Step& operator=(const Step&) = delete;
    Step& operator=(Step&&) noexcept = delete;
    ~Step() override = default;

    virtual int initialize();

    virtual int analyze() = 0;

    void set_domain(const weak_ptr<DomainBase>&);
    [[nodiscard]] const weak_ptr<DomainBase>& get_domain() const;

    void set_factory(const shared_ptr<Factory<double>>&);
    [[nodiscard]] const shared_ptr<Factory<double>>& get_factory() const;

    void set_solver_tag(unsigned);
    void set_solver(const shared_ptr<Solver>&);
    [[nodiscard]] const shared_ptr<Solver>& get_solver() const;

    void set_converger_tag(unsigned);
    void set_converger(const shared_ptr<Converger>&);
    [[nodiscard]] const shared_ptr<Converger>& get_converger() const;

    void set_integrator_tag(unsigned);
    void set_integrator(const shared_ptr<Integrator>&);
    [[nodiscard]] const shared_ptr<Integrator>& get_integrator() const;

    void set_time_period(double);
    void set_time_left(double);
    [[nodiscard]] double get_time_period() const;
    [[nodiscard]] double get_time_left() const;

    void set_ini_step_size(double);
    void set_min_step_size(double);
    void set_max_step_size(double);
    void set_max_substep(unsigned);
    void set_system_solver(SolverType);
    void set_sub_system_solver(SolverType);
    void set_precision(Precision);
    void set_tolerance(double);
    void set_refinement(std::uint8_t);
    void set_lis_option(std::string_view);
#ifdef SUANPAN_MAGMA
    void set_magma_option(const magma_dopts& magma_opt) { magma_setting = magma_opt; }
#endif

    [[nodiscard]] double get_ini_step_size() const;
    [[nodiscard]] double get_min_step_size() const;
    [[nodiscard]] double get_max_step_size() const;
    [[nodiscard]] unsigned get_max_substep() const;

    [[nodiscard]] bool is_fixed_step_size() const;
    void set_fixed_step_size(bool);

    [[nodiscard]] bool is_symm() const;
    [[nodiscard]] bool is_band() const;
    [[nodiscard]] bool is_sparse() const;
    void set_symm(bool) const;
    void set_band(bool) const;
    void set_sparse(bool) const;
};

#endif

//! @}
