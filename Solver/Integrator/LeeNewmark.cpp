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

#include "LeeNewmark.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

uword LeeNewmark::get_total_size() const { return n_damping * n_block + n_block; }

void LeeNewmark::update_stiffness() const {
	if(StorageScheme::SPARSE == factory->get_storage_scheme())
		for(uword I = 0, J = n_block; I < n_damping; ++I, J += n_block) {
			auto row = current_mass->triplet_mat.row_mem();
			auto col = current_mass->triplet_mat.col_mem();
			auto val = current_mass->triplet_mat.val_mem();
			for(size_t O = 0; O < current_mass->triplet_mat.c_size; ++O) {
				const auto K = row[O], L = col[O], M = K + J, N = L + J;
				stiffness->at(M, L) = C1 * (stiffness->at(K, N) = -(stiffness->at(M, N) = mass_coef(I) * val[O]));
			}
			row = current_stiffness->triplet_mat.row_mem();
			col = current_stiffness->triplet_mat.col_mem();
			val = current_stiffness->triplet_mat.val_mem();
			for(size_t K = 0; K < current_stiffness->triplet_mat.c_size; ++K) stiffness->at(row[K] + J, col[K] + J) = stiffness_coef(I) * val[K];
		}
	else
		for(uword I = 0, J = n_block; I < n_damping; ++I, J += n_block)
			for(unsigned K = 0; K < n_block; ++K) {
				const auto M = K + J;
				for(unsigned L = 0; L < n_block; ++L) {
					const auto N = L + J;
					auto t_val = current_mass->operator()(K, L);
					if(0. != t_val) stiffness->at(M, L) = C1 * (stiffness->at(K, N) = -(stiffness->at(M, N) = mass_coef(I) * t_val));
					t_val = current_stiffness->operator()(K, L);
					if(0. != t_val) stiffness->at(M, N) = stiffness_coef(I) * t_val;
				}
			}
}

void LeeNewmark::update_residual() const {
	const auto& t_vel = factory->get_trial_velocity();

	auto& t_residual = access::rw(residual);

	for(uword I = 0, J = n_block, K = J + n_block - 1; I < n_damping; ++I, J += n_block, K += n_block) {
		const vec n_internal(access::rwp(&trial_internal(J)), n_block, false, true);
		t_residual.rows(J, K) = current_mass * vec(t_vel - n_internal) * mass_coef(I) - current_stiffness * n_internal * stiffness_coef(I);
	}
}

LeeNewmark::LeeNewmark(const unsigned T, vec&& X, vec&& F, const double A, const double B)
	: LeeNewmarkBase(T, A, B)
	, mass_coef(4. * X % F)
	, stiffness_coef(4. * X / F)
	, CM(4. * dot(X, F)) {}

int LeeNewmark::initialize() {
	if(SUANPAN_SUCCESS != LeeNewmarkBase::initialize()) return SUANPAN_FAIL;

	current_mass = factory->get_mass()->make_copy();
	current_stiffness = factory->get_stiffness()->make_copy();

	return SUANPAN_SUCCESS;
}

int LeeNewmark::process_constraint() {
	const auto& D = get_domain().lock();

	// process constraint for the first time to obtain proper stiffness
	if(SUANPAN_SUCCESS != Integrator::process_constraint()) return SUANPAN_FAIL;

	auto& t_stiff = get_stiffness(factory);
	auto& t_mass = get_mass(factory);

	t_stiff->csc_condense();

	if(first_iteration) {
		stiffness->triplet_mat.resize((4 * n_damping + 2) * t_stiff->n_elem);
		stiffness->zeros();

		// current_mass.swap(t_mass);
		// D->assemble_current_mass();
		// current_mass.swap(t_mass);

		current_mass = t_mass->make_copy();
	}
	else {
		// if not first iteration
		// erase the tangent stiffness entries
		if(!stiffness->triplet_mat.csc_sort()) return SUANPAN_FAIL;

		const auto& row = stiffness->triplet_mat.row_mem();
		const auto& col = stiffness->triplet_mat.col_mem();
		const auto& val = stiffness->triplet_mat.val_idx;

		for(size_t I = 0; I < stiffness->triplet_mat.c_size; ++I) {
			// quit if current column is beyond the original size of matrix
			if(col[I] >= n_block) break;
			// erase existing entries if fall in intact stiffness matrix
			if(row[I] < n_block) val[I] = 0.;
		}
	}

	t_stiff += C1 * CM * current_mass;

	t_stiff->csc_condense();

	// check in tangent stiffness
	if(StorageScheme::SPARSE == factory->get_storage_scheme()) stiffness->triplet_mat += t_stiff->triplet_mat;
	else for(unsigned I = 0; I < n_block; ++I) for(unsigned J = 0; J < n_block; ++J) if(const auto t_val = t_stiff->operator()(I, J); 0. != t_val) stiffness->at(I, J) = t_val;

	if(first_iteration) {
		// for the first iteration of each substep
		// store current stiffness to be used in the whole substep
		// check in constant terms that does not change in the substep
		current_stiffness.swap(t_stiff);
		D->assemble_current_stiffness();
		if(SUANPAN_SUCCESS != Integrator::process_constraint()) return SUANPAN_FAIL;
		t_stiff->csc_condense();
		current_stiffness.swap(t_stiff);

		update_stiffness();

		first_iteration = false;
	}

	update_residual();

	return SUANPAN_SUCCESS;
}

void LeeNewmark::assemble_resistance() {
	const auto& D = get_domain().lock();
	const auto& W = factory;

	D->assemble_resistance();
	D->assemble_inertial_force();
	// consider independent viscous device
	D->assemble_damping_force();

	if(nullptr != current_mass) {
		vec internal_velocity = CM * W->get_trial_velocity();
		for(uword I = 0, J = n_block; I < n_damping; ++I, J += n_block) {
			const vec n_internal(&trial_internal(J), n_block, false, true);
			internal_velocity -= mass_coef(I) * n_internal;
		}
		W->update_trial_damping_force(W->get_trial_damping_force() + current_mass * internal_velocity);
	}

	W->set_sushi(W->get_trial_resistance() + W->get_trial_damping_force() + W->get_trial_inertial_force());
}

void LeeNewmark::print() {
	suanpan_info("A Newmark solver using Lee's damping model. DOI: 10.1016/j.jsv.2020.115312\n");
	const vec X = .25 * sqrt(mass_coef % stiffness_coef);
	const vec F = sqrt(mass_coef / stiffness_coef);
	for(auto I = 0llu; I < n_damping; ++I) suanpan_info("\tdamping ratio: %.4f\tfrequency: %.4f\n", X(I), F(I));
}
