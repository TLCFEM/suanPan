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

	new_feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
	fpm[0] = 1;
#endif

	std::vector output(4, 0);
	std::vector input(4, 0.);
	input[1] = 0.;     // centre
	input[2] = radius; // radius

	output[1] = eigen_num;

	int N = static_cast<int>(W->get_size());

	auto M = static_cast<int>(eigen_num);
	std::vector R(M, 0.);
	std::vector E(M, 0.);
	M *= N;
	std::vector X(M, 0.);

	char UPLO = 'F';

	if(const auto scheme = W->get_storage_scheme(); StorageScheme::FULL == scheme) new_dfeast_sygv_(&UPLO, &N, stiffness->memptr(), &N, mass->memptr(), &N, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
	else if(StorageScheme::SPARSE == scheme) {
		auto fs = std::async([&]() { return csr_form<double, int>(stiffness->triplet_mat, 1); });
		auto fm = std::async([&]() { return csr_form<double, int>(mass->triplet_mat, 1); });

		const auto t_stiff = fs.get();
		const auto t_mass = fm.get();

		new_dfeast_scsrgv_(&UPLO, &N, t_stiff.val_idx, t_stiff.row_ptr, t_stiff.col_idx, t_mass.val_idx, t_mass.row_ptr, t_mass.col_idx, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
	}
	else if(StorageScheme::BAND == scheme || StorageScheme::BANDSYMM == scheme) {
		fpm[41] = 0;

		unsigned l, u;
		W->get_bandwidth(l, u);
		auto KL = static_cast<int>(l);
		const auto KU = static_cast<int>(u);
		auto LD = KL + KU + 1;

		new_dfeast_sbgv_(&UPLO, &N, &KL, stiffness->memptr(), &LD, &KL, mass->memptr(), &LD, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
	}
	else throw;

	if(0 != output[3]) {
		suanpan_error("error code %d recieved from FEAST solver.\n", output[3]);
		return SUANPAN_FAIL;
	}

	auto& eigval = get_eigenvalue(W);
	eigval.set_size(output[2]);

#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, eigval.n_elem, [&](const uword I) { eigval(I) = E[I]; });
#else
	for(uword I = 0; I < eigval.n_elem; ++I) eigval(I) = E[I];
#endif

	auto& eigvec = get_eigenvector(W);
	eigvec.resize(N, output[2]);

#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, eigvec.n_elem, [&](const uword I) { eigvec(I) = X[I]; });
#else
	for(uword I = 0; I < eigvec.n_elem; ++I) eigvec(I) = X[I];
#endif

	return SUANPAN_SUCCESS;
}

int FEAST::quadratic_solve(const shared_ptr<LongFactory>& W) const {
	std::vector fpm(64, 0);

	new_feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
	fpm[0] = 1;
#endif

	int N = static_cast<int>(W->get_size());

	std::vector output(4, 0);
	std::vector input(4, 0.);
	input[0] = radius; // centre
	input[1] = 0.;     // centre
	input[2] = radius; // radius

	output[1] = eigen_num;

	auto M = 2 * static_cast<int>(eigen_num);
	std::vector R(M, 0.);
	std::vector E(M, 0.);
	M *= 2 * N;
	std::vector X(M, 0.);

	int P = 2;

	const auto& mass = W->get_mass();
	const auto& damping = W->get_damping();
	const auto& stiffness = W->get_stiffness();

	auto& c_stiff = stiffness->triplet_mat;
	auto& c_damping = damping->triplet_mat;
	auto& c_mass = mass->triplet_mat;

	if(StorageScheme::SPARSE != W->get_storage_scheme()) {
		const auto NN = N * N;

		auto fk = std::async([&]() {
			c_stiff.init(N, N, NN);
			for(auto I = 0; I < N; ++I) for(auto J = 0; J < N; ++J) c_stiff.at(I, J) = stiffness->operator()(I, J);
		});
		auto fd = std::async([&]() {
			c_damping.init(N, N, NN);
			for(auto I = 0; I < N; ++I) for(auto J = 0; J < N; ++J) c_damping.at(I, J) = damping->operator()(I, J);
		});
		auto fm = std::async([&]() {
			c_mass.init(N, N, NN);
			for(auto I = 0; I < N; ++I) for(auto J = 0; J < N; ++J) c_mass.at(I, J) = mass->operator()(I, J);
		});

		fk.get();
		fd.get();
		fm.get();
	}

	auto fk = std::async([&]() { return csr_form<double, int>(c_stiff, 1); });
	auto fd = std::async([&]() { return csr_form<double, int>(c_damping, 1); });
	auto fm = std::async([&]() { return csr_form<double, int>(c_mass, 1); });

	const auto t_stiff = fk.get();
	const auto t_damping = fd.get();
	const auto t_mass = fm.get();

	size_t n_elem1 = std::max(std::max(t_stiff.c_size, t_damping.c_size), t_mass.c_size);
	size_t n_elem2 = n_elem1 + n_elem1;
	size_t n_elem3 = n_elem1 + n_elem2;
	size_t n_size1 = N;
	size_t n_size2 = n_size1 + n_size1;
	size_t n_size3 = n_size2 + n_size1;

	std::vector A(n_elem3, 0.);
	std::vector JA(n_elem3, 0);
	std::vector IA(n_size3, 0);

#ifdef SUANPAN_MT
	tbb::parallel_for(0, t_stiff.c_size, [&](const int I) {
		A[I] = t_stiff.val_idx[I];
		JA[I] = t_stiff.col_idx[I];
	});
	tbb::parallel_for(0, t_damping.c_size, [&](const int I) {
		A[I + n_elem1] = t_damping.val_idx[I];
		JA[I + n_elem1] = t_damping.col_idx[I];
	});
	tbb::parallel_for(0, t_mass.c_size, [&](const int I) {
		A[I + n_elem2] = t_mass.val_idx[I];
		JA[I + n_elem2] = t_mass.col_idx[I];
	});

	tbb::parallel_for(0, N, [&](const int I) {
		JA[I] = t_stiff.row_ptr[I];
		JA[I + n_size1] = t_damping.row_ptr[I];
		JA[I + n_size2] = t_mass.row_ptr[I];
	});
#else
	std::copy_n(t_stiff.val_idx, t_stiff.c_size, A.data());
	std::copy_n(t_damping.val_idx, t_damping.c_size, A.data() + n_elem1);
	std::copy_n(t_mass.val_idx, t_mass.c_size, A.data() + n_elem2);

	std::copy_n(t_stiff.col_idx, t_stiff.c_size, JA.data());
	std::copy_n(t_damping.col_idx, t_damping.c_size, JA.data() + n_elem1);
	std::copy_n(t_mass.col_idx, t_mass.c_size, JA.data() + n_elem2);

	std::copy_n(t_stiff.row_ptr, N, IA.data());
	std::copy_n(t_damping.row_ptr, N, IA.data() + n_size1);
	std::copy_n(t_mass.row_ptr, N, IA.data() + n_size2);
#endif

	new_dfeast_gcsrpev_(&P, &N, A.data(), IA.data(), JA.data(), fpm.data(), &input[3], &output[0], &input[0], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);

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
