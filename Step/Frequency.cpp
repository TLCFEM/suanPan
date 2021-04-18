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

#include "Frequency.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Arnoldi.h>

Frequency::Frequency(const unsigned T, const unsigned N)
	: Step(T, 0.)
	, eigen_number(N) {}

int Frequency::initialize() {
	if(sparse_mat) factory->set_storage_scheme(StorageScheme::SPARSE);
	else if(symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BANDSYMM);
	else if(!symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BAND);
	else if(symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::SYMMPACK);
	else if(!symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::FULL);

	factory->set_analysis_type(AnalysisType::EIGEN);
	factory->set_precision(precision);
	factory->set_tolerance(tolerance);
	factory->set_solver(system_solver);

	const auto& t_domain = database.lock();

	modifier = make_shared<Integrator>();
	if(nullptr == solver) solver = make_shared<Arnoldi>(0, eigen_number);

	modifier->set_domain(t_domain);
	solver->set_integrator(modifier);

	return t_domain->restart();
}

int Frequency::analyze() {
	auto& G = get_integrator();

	if(SUANPAN_SUCCESS != solver->analyze()) {
		suanpan_warning("fail to decompose the system, try to increase the number of eigen values.\n");
		return SUANPAN_SUCCESS;
	}

	G->record();

	return SUANPAN_SUCCESS;
}

void Frequency::set_eigen_number(const unsigned N) const { access::rw(eigen_number) = N; }

unsigned Frequency::get_eigen_number() const { return eigen_number; }
