/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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

#include "Buckle.h"
#include <Converger/AbsIncreDisp.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Newton.h>
#include <Toolbox/arpack_wrapper.h>

Buckle::Buckle(const unsigned T)
	: Static(T, 1.) {}

int Buckle::initialize() {
	const auto& t_domain = database.lock();

	modifier = make_shared<Integrator>();

	if(nullptr == tester) tester = make_shared<AbsIncreDisp>();
	if(nullptr == solver) solver = make_shared<Newton>();

	tester->set_domain(t_domain);
	modifier->set_domain(t_domain);
	solver->set_converger(tester);
	solver->set_integrator(modifier);

	if(sparse_mat) factory->set_storage_scheme(StorageScheme::SPARSE);
	else if(symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BANDSYMM);
	else if(!symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BAND);
	else if(symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::SYMMPACK);
	else if(!symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::FULL);

	factory->set_analysis_type(AnalysisType::BUCKLE);
	factory->set_precision(precision);
	factory->set_tolerance(tolerance);
	factory->set_solver(system_solver);

	return t_domain->restart();
}

int Buckle::analyze() {
	if(Static::analyze() == SUANPAN_FAIL) return SUANPAN_FAIL;

	const auto& D = get_domain().lock();
	auto& G = get_integrator();
	auto& W = get_factory();

	// assemble stiffness and geometry as they may be modified in solver
	D->assemble_trial_stiffness();
	D->assemble_trial_geometry();

	// now need to apply constraints on both stiffness and geometry
	if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;
	// swap stiffness and geometry
	access::rw(W->get_stiffness()).swap(access::rw(W->get_geometry()));
	// apply constraints again on geometry part
	if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;
	// swap back
	access::rw(W->get_stiffness()).swap(access::rw(W->get_geometry()));

	if(eig_solve(get_eigenvalue(W), get_eigenvector(W), W->get_stiffness(), W->get_geometry()) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	suanpan_info("\nbuckling load multiplier: %.8E.\n", W->get_eigenvalue().at(0));

	// record response
	G->record();

	return SUANPAN_SUCCESS;
}
