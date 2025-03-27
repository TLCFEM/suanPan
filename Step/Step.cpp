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

#include "Step.h"
#include <Converger/RelIncreDisp.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Solver.h>

void Step::configure_storage_scheme() const {
    if(sparse_mat) factory->set_storage_scheme(symm_mat ? StorageScheme::SPARSESYMM : StorageScheme::SPARSE);
    else if(symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BANDSYMM);
    else if(!symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BAND);
    else if(symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::SYMMPACK);
    else if(!symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::FULL);
}

Step::Step(const unsigned T, const double P)
    : Tag(T)
    , time_period(P) {}

int Step::initialize() {
    const auto t_domain = database.lock();

    if(sparse_mat) {
        // LAPACK and SPIKE are for dense only
        if(SolverType::LAPACK == system_solver || SolverType::SPIKE == system_solver) system_solver = SolverType::SUPERLU;
    }
    else if(!symm_mat && band_mat) {
        // only LAPACK and SPIKE are supported
        if(SolverType::LAPACK != system_solver && SolverType::SPIKE != system_solver) system_solver = SolverType::LAPACK;
    }
    else if(!symm_mat && !band_mat) {
        // only LAPACK and CUDA solvers are supported
        if(SolverType::LAPACK != system_solver && SolverType::CUDA != system_solver) system_solver = SolverType::LAPACK;
    }

    if(converger_tag != 0 && t_domain->find_converger(converger_tag)) tester = t_domain->get_converger(converger_tag);
    else if(const auto [t_tag, step_tag] = t_domain->get_current_converger_tag(); t_tag != 0 && step_tag <= get_tag()) tester = t_domain->get_current_converger();

    if(integrator_tag != 0 && t_domain->find_integrator(integrator_tag)) modifier = t_domain->get_integrator(integrator_tag);
    else if(const auto [t_tag, step_tag] = t_domain->get_current_integrator_tag(); t_tag != 0 && step_tag <= get_tag()) modifier = t_domain->get_current_integrator();

    if(solver_tag != 0 && t_domain->find_solver(solver_tag)) solver = t_domain->get_solver(solver_tag);
    else if(const auto [t_tag, step_tag] = t_domain->get_current_solver_tag(); t_tag != 0 && step_tag <= get_tag()) solver = t_domain->get_current_solver();

    if(tester == nullptr) tester = make_shared<RelIncreDisp>();

    factory = t_domain->get_factory();

    factory->set_solver_type(system_solver);
    factory->set_solver_setting(system_setting);
    factory->set_sub_solver_type(sub_system_solver);
#ifdef SUANPAN_MAGMA
    factory->set_solver_setting(magma_setting);
#endif

    return 0;
}

void Step::set_domain(const weak_ptr<DomainBase>& D) { database = D; }

const weak_ptr<DomainBase>& Step::get_domain() const { return database; }

void Step::set_factory(const shared_ptr<Factory<double>>& F) { factory = F; }

const shared_ptr<Factory<double>>& Step::get_factory() const { return factory; }

void Step::set_solver_tag(const unsigned T) { solver_tag = T; }

void Step::set_solver(const shared_ptr<Solver>& S) { solver = S; }

const shared_ptr<Solver>& Step::get_solver() const { return solver; }

void Step::set_converger_tag(const unsigned T) { converger_tag = T; }

void Step::set_converger(const shared_ptr<Converger>& C) { tester = C; }

const shared_ptr<Converger>& Step::get_converger() const { return tester; }

void Step::set_integrator_tag(const unsigned T) { integrator_tag = T; }

void Step::set_integrator(const shared_ptr<Integrator>& G) { modifier = G; }

const shared_ptr<Integrator>& Step::get_integrator() const { return modifier; }

void Step::set_time_period(const double T) {
    if(suanpan::approx_equal(time_period, T)) return;
    time_left = time_period = T;
    const auto t_iteration = static_cast<int>(floor(time_period / ini_step_size)) + 1;
    if(t_iteration <= static_cast<int>(max_substep) || 0 == max_substep) return;
    if(t_iteration > static_cast<int>(std::numeric_limits<unsigned>::max())) {
        suanpan_warning("Too many iterations.\n");
        set_max_substep(std::numeric_limits<unsigned>::max());
    }
    else set_max_substep(t_iteration);
}

void Step::set_time_left(const double T) { time_left = T; }

double Step::get_time_period() const { return time_period; }

double Step::get_time_left() const { return time_left; }

void Step::set_ini_step_size(const double T) {
    if(suanpan::approx_equal(ini_step_size, T)) return;
    if(time_period > 0.) {
        // for time control
        ini_step_size = T > time_period ? time_period : T;
        if(const auto t_iteration = static_cast<int>(floor(time_period / ini_step_size)) + 1; t_iteration > static_cast<int>(max_substep) && max_substep != 0) set_max_substep(t_iteration);
    }
    else ini_step_size = fabs(T); // for arc-length control
}

void Step::set_min_step_size(const double T) { min_step_size = T; }

void Step::set_max_step_size(const double T) { max_step_size = T; }

void Step::set_max_substep(const unsigned M) { max_substep = M; }

void Step::set_system_solver(const SolverType P) { system_solver = P; }

void Step::set_sub_system_solver(const SolverType P) { sub_system_solver = P; }

void Step::set_precision(const Precision P) { system_setting.precision = P; }

void Step::set_tolerance(const double T) { system_setting.tolerance = T; }

void Step::set_refinement(const std::uint8_t T) { system_setting.iterative_refinement = T; }

void Step::set_lis_option(const std::string_view T) { system_setting.lis_options = T; }

double Step::get_ini_step_size() const { return ini_step_size; }

double Step::get_min_step_size() const { return min_step_size; }

double Step::get_max_step_size() const { return max_step_size; }

unsigned Step::get_max_substep() const { return max_substep; }

bool Step::is_fixed_step_size() const { return fixed_step_size; }

void Step::set_fixed_step_size(const bool B) { fixed_step_size = B; }

bool Step::is_symm() const { return symm_mat; }

bool Step::is_band() const { return band_mat; }

bool Step::is_sparse() const { return sparse_mat; }

void Step::set_symm(const bool B) const { access::rw(symm_mat) = B; }

void Step::set_band(const bool B) const { access::rw(band_mat) = B; }

void Step::set_sparse(const bool B) const { access::rw(sparse_mat) = B; }
