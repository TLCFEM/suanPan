////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "FEAST.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Toolbox/arpack_wrapper.h>

#ifndef SUANPAN_MKL
#include <feast/feast.h>
#endif

FEAST::FEAST(const unsigned T, const unsigned N, const double R)
	: Solver(T)
	, eigen_num(N)
	, radius(R) {}

int FEAST::analyze() {
	auto& G = get_integrator();
	const auto& D = G->get_domain().lock();
	auto& W = D->get_factory();

	D->assemble_trial_mass();
	D->assemble_trial_stiffness();

	if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
	if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;

	const auto& stiffness = W->get_stiffness();
	const auto& mass = W->get_mass();

	const auto scheme = W->get_storage_scheme();

	if(StorageScheme::SYMMPACK == scheme) return eig_solve(get_eigenvalue(W), get_eigenvector(W), stiffness, mass, eigen_num, "SM");

#ifdef SUANPAN_MKL
	return eig_solve(get_eigenvalue(W), get_eigenvector(W), stiffness, mass, eigen_num, "SM");
#else
	podarray<int> fpm(64);

	feastinit_(fpm.mem);

	fpm(0) = 1;
	fpm(14) = 1;

	if(StorageScheme::FULL == scheme) {
		podarray<int> output(4);
		podarray<double> input(4);
		input(0) = radius; // centre
		input(1) = 0.;     // centre
		input(2) = radius; // radius

		int M = 2 * static_cast<int>(eigen_num);
		const podarray<double> R(M);
		const podarray<double> E(M);
		int N = static_cast<int>(W->get_size());
		M *= N;
		const podarray<double> X(M);

		output(1) = eigen_num;

		dfeast_gegv_(&N, stiffness->memptr(), &N, mass->memptr(), &N, fpm.mem, &input(3), &output(0), input.mem, &input(2), &output(1), E.mem, X.mem, &output(2), R.mem, &output(3));

		if(0 != output(3)) {
			suanpan_error("error code %d recieved from FEAST solver.\n", output(3));
			return SUANPAN_FAIL;
		}

		auto& eigval = get_eigenvalue(W);
		eigval.set_size(output(2));

		for(uword I = 0; I < eigval.n_elem; ++I) eigval(I) = E(2 * I);

		auto& eigvec = get_eigenvector(W);
		eigvec.resize(N, output(2));

		for(uword I = 0; I < eigvec.n_elem; ++I) eigvec(I) = X(2 * I);
	}
#endif

	return SUANPAN_SUCCESS;
}

void FEAST::print() { suanpan_info("An eigen solver using FEAST solver.\n"); }
