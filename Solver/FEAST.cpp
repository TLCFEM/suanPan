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
#include <feast/feast.h>

int FEAST::linear_solve(const shared_ptr<LongFactory>& W) const {
	const auto& mass = W->get_mass();
	const auto& stiffness = W->get_stiffness();

	std::vector fpm(64, 0);

	feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
	fpm[0] = 1;
#endif

	int N = static_cast<int>(W->get_size());

	std::vector output(4, 0);
	std::vector input(4, 0.);
	input[1] = 0.;     // centre
	input[2] = radius; // radius

	auto M = static_cast<int>(eigen_num);
	std::vector R(M, 0.);
	std::vector E(M, 0.);
	M *= N;
	std::vector X(M, 0.);

	output[1] = eigen_num;

	char UPLO = 'F';

	if(const auto scheme = W->get_storage_scheme(); StorageScheme::FULL == scheme) dfeast_sygv_(&UPLO, &N, stiffness->memptr(), &N, mass->memptr(), &N, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
	else if(StorageScheme::SPARSE == scheme) {
		const csr_form<double, int> t_stiff(stiffness->triplet_mat, 1);
		const csr_form<double, int> t_mass(mass->triplet_mat, 1);

		dfeast_scsrgv_(&UPLO, &N, t_stiff.val_idx, t_stiff.row_ptr, t_stiff.col_idx, t_mass.val_idx, t_mass.row_ptr, t_mass.col_idx, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
	}
	else if(StorageScheme::BAND == scheme || StorageScheme::BANDSYMM == scheme) {
		fpm[41] = 0;

		unsigned l, u;
		W->get_bandwidth(l, u);
		auto KL = static_cast<int>(l);
		const auto KU = static_cast<int>(u);
		auto LD = KL + KU + 1;

		dfeast_sbgv_(&UPLO, &N, &KL, stiffness->memptr(), &LD, &KL, mass->memptr(), &LD, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
	}
	else throw;

	if(0 != output[3]) {
		suanpan_error("error code %d recieved from FEAST solver.\n", output[3]);
		return SUANPAN_FAIL;
	}

	auto& eigval = get_eigenvalue(W);
	eigval.set_size(output[2]);

	for(uword I = 0; I < eigval.n_elem; ++I) eigval(I) = E[I];

	auto& eigvec = get_eigenvector(W);
	eigvec.resize(N, output[2]);

	for(uword I = 0; I < eigvec.n_elem; ++I) eigvec(I) = X[I];

	return SUANPAN_SUCCESS;
}

int FEAST::quadratic_solve(const shared_ptr<LongFactory>& W) const {
	std::vector fpm(64, 0);

	feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
	fpm[0] = 1;
#endif

	int N = static_cast<int>(W->get_size());

	std::vector output(4, 0);
	std::vector input(4, 0.);
	input[0] = radius; // centre
	input[1] = 0.;     // centre
	input[2] = radius; // radius

	auto M = 2 * static_cast<int>(eigen_num);
	std::vector R(M, 0.);
	std::vector E(M, 0.);
	M *= 2 * N;
	std::vector X(M, 0.);

	output[1] = eigen_num;

	int P = 2;

	const auto& mass = W->get_mass();
	const auto& damping = W->get_damping();
	const auto& stiffness = W->get_stiffness();

	auto& c_stiff = stiffness->triplet_mat;
	auto& c_damping = damping->triplet_mat;
	auto& c_mass = mass->triplet_mat;

	if(StorageScheme::SPARSE != W->get_storage_scheme()) {
		const auto NN = N * N;
		c_stiff.init(N, N, NN);
		c_damping.init(N, N, NN);
		c_mass.init(N, N, NN);

		for(auto I = 0; I < N; ++I)
			for(auto J = 0; J < N; ++J) {
				c_stiff.at(I, J) = stiffness->operator()(I, J);
				c_damping.at(I, J) = damping->operator()(I, J);
				c_mass.at(I, J) = mass->operator()(I, J);
			}
	}

	const csr_form<double, int> t_stiff(c_stiff, 1);
	const csr_form<double, int> t_damping(c_damping, 1);
	const csr_form<double, int> t_mass(c_mass, 1);

	auto n_elem = std::max(std::max(t_stiff.c_size, t_damping.c_size), t_mass.c_size);

	std::vector A(3llu * n_elem, 0.);
	std::vector JA(3llu * n_elem, 0);
	std::vector IA(3llu * N, 0);

	std::copy_n(t_stiff.val_idx, t_stiff.c_size, A.data());
	std::copy_n(t_damping.val_idx, t_damping.c_size, A.data() + n_elem);
	std::copy_n(t_mass.val_idx, t_mass.c_size, A.data() + 2llu * n_elem);

	std::copy_n(t_stiff.col_idx, t_stiff.c_size, JA.data());
	std::copy_n(t_damping.col_idx, t_damping.c_size, JA.data() + n_elem);
	std::copy_n(t_mass.col_idx, t_mass.c_size, JA.data() + 2llu * n_elem);

	std::copy_n(t_stiff.row_ptr, N, IA.data());
	std::copy_n(t_damping.row_ptr, N, IA.data() + N);
	std::copy_n(t_mass.row_ptr, N, IA.data() + 2llu * N);

	dfeast_gcsrpev_(&P, &N, A.data(), IA.data(), JA.data(), fpm.data(), &input[3], &output[0], &input[0], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);

	if(0 != output[3]) {
		suanpan_error("error code %d recieved from FEAST solver.\n", output[3]);
		return SUANPAN_FAIL;
	}

	auto& eigval = get_eigenvalue(W);
	eigval.set_size(output[2]);

	for(uword I = 0; I < eigval.n_elem; ++I) eigval(I) = E[2 * I];

	auto& eigvec = get_eigenvector(W);
	eigvec.resize(N, output[2]);

	for(uword I = 0; I < eigvec.n_elem; ++I) eigvec(I) = X[2 * I];

	return SUANPAN_SUCCESS;
}

FEAST::FEAST(const unsigned T, const unsigned N, const double R, const bool Q)
	: Solver(T)
	, quadratic(Q)
	, eigen_num(N)
	, radius(R) {}

int FEAST::initialize() {
	if(SUANPAN_SUCCESS != Solver::initialize()) return SUANPAN_FAIL;

	auto& G = get_integrator();
	const auto& D = G->get_domain().lock();
	auto& W = D->get_factory();

	if(const auto scheme = W->get_storage_scheme(); StorageScheme::SYMMPACK == scheme) {
		suanpan_error("FEAST solver does not support symmetric pack storage.\n");

		return SUANPAN_FAIL;
	}
	else if((StorageScheme::BAND == scheme || StorageScheme::BANDSYMM == scheme) && SolverType::SPIKE != W->get_solver()) {
		suanpan_error("SPIKE system solver needs to be used for banded storage.\n");

		return SUANPAN_FAIL;
	}

	return SUANPAN_SUCCESS;
}

int FEAST::analyze() {
	auto& G = get_integrator();
	const auto& D = G->get_domain().lock();
	auto& W = D->get_factory();

	D->assemble_trial_mass();
	D->assemble_trial_stiffness();
	if(quadratic) D->assemble_trial_damping();

	if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
	if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;

	return quadratic ? quadratic_solve(W) : linear_solve(W);
}

void FEAST::print() { suanpan_info("An eigen solver using FEAST solver.\n"); }
