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

#include "LeeNewmarkFull.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

/**
 * \brief compute an approximation of the number of additional blocks
 * \return amplifier
 */
uword LeeNewmarkFull::get_amplifier() const {
	// TODO: CHECK ACCURACY

	auto n_size = 2llu;

	for(const auto& [t, p, zeta, omega] : damping_mode)
		if(Type::T0 == t) n_size += 5llu;
		else if(Type::T1 == t) n_size += 5llu + 6llu * static_cast<uword>(p.front());
		else if(Type::T2 == t) n_size += 4llu + 5llu * static_cast<uword>(.5 * (p(0) + p(1) - 1.));
		else if(Type::T3 == t) n_size += 9llu;

	return n_size;
}

/**
 * \brief compute size of final global effective stiffness
 * \return the size of final global effective stiffness
 */
uword LeeNewmarkFull::get_total_size() const {
	auto n_size = 1llu;

	for(const auto& [t, p, zeta, omega] : damping_mode)
		if(Type::T0 == t) n_size += 1llu;
		else if(Type::T1 == t) n_size += 2llu * static_cast<uword>(p.front()) + 1llu;
		else if(Type::T2 == t) n_size += static_cast<uword>(p(0) + p(1)) + 1llu;
		else if(Type::T3 == t) n_size += 2llu;

	return n_size * n_block;
}

void LeeNewmarkFull::update_stiffness() const {
	auto IDX = n_block;

	// ! make sure global stiffness only holds unrolled damping matrix when exit
	stiffness->zeros();

	for(const auto& [t, p, zeta, omega] : damping_mode)
		if(const auto mass_coef = 4. * zeta * omega * C1, stiffness_coef = 4. * zeta / omega * C1; Type::T0 == t) assemble_by_mode_zero(IDX, mass_coef, stiffness_coef);
		else if(Type::T1 == t) assemble_by_mode_one(IDX, mass_coef, stiffness_coef, p.front());
		else if(Type::T2 == t) assemble_by_mode_two(IDX, mass_coef, stiffness_coef, p(0), p(1));
		else if(Type::T3 == t) assemble_by_mode_three(IDX, mass_coef, stiffness_coef, p.front());

	stiffness->csc_condense();
}

void LeeNewmarkFull::update_residual() const {
	// ! please make sure global stiffness holds unrolled damping matrix only when calling this method

	// ! only account for residual due to damping matrix
	// ! will be completed after calling the method to get residual
	vec trial_vel = -trial_internal;
	trial_vel.head(n_block) = factory->get_trial_velocity() / -C1;
	access::rw(residual) = stiffness * trial_vel;
	// ? may potentially improve performance
	// access::rw(residual) = csr_form<double>(stiffness->triplet_mat) * trial_vel;

	// ! check in damping force
	auto fa = std::async([&]() { get_trial_damping_force(factory) -= residual.head(n_block); });
	auto fb = std::async([&]() { get_incre_damping_force(factory) -= residual.head(n_block); });
	// ! update left hand side
	auto fc = std::async([&]() { get_sushi(factory) -= residual.head(n_block); });

	fa.get();
	fb.get();
	fc.get();
}

void LeeNewmarkFull::assemble_by_mode_zero(uword& current_pos, const double mass_coef, const double stiffness_coef) const {
	const auto I = current_pos;

	auto row = current_mass.row_mem();
	auto col = current_mass.col_mem();
	auto val = current_mass.val_mem();
	for(index_tm J = 0; J < current_mass.c_size; ++J) {
		const auto K = row[J], L = col[J], M = K + I, N = L + I;
		stiffness->at(K, N) = stiffness->at(M, L) = -(stiffness->at(K, L) = stiffness->at(M, N) = mass_coef * val[J]);
	}
	row = current_stiffness.row_mem();
	col = current_stiffness.col_mem();
	val = current_stiffness.val_mem();
	for(index_ts J = 0; J < current_stiffness.c_size; ++J) stiffness->at(row[J] + I, col[J] + I) = stiffness_coef * val[J];

	current_pos += n_block;
}

void LeeNewmarkFull::assemble_by_mode_one(uword& current_pos, const double mass_coef, const double stiffness_coef, double order) const {
	const auto mass_coefs = .5 * mass_coef;           // eq. 10
	const auto stiffness_coefs = .5 * stiffness_coef; // eq. 10

	const auto& m_row = current_mass.row_mem();
	const auto& m_col = current_mass.col_mem();
	const auto& m_val = current_mass.val_mem();
	const auto& s_row = current_stiffness.row_mem();
	const auto& s_col = current_stiffness.col_mem();
	const auto& s_val = current_stiffness.val_mem();

	auto I = 0llu;
	auto J = current_pos;
	auto K = current_pos += n_block;
	auto L = current_pos += n_block;
	auto M = current_pos += n_block;

	while(order > 1.5) {
		// eq. 61 
		for(index_tm N = 0; N < current_mass.c_size; ++N) {
			const auto O = m_row[N], P = m_col[N], Q = O + J, R = P + J;
			stiffness->at(O + I, R) = stiffness->at(Q, P + I) = mass_coef * m_val[N];
			stiffness->at(Q, P + K) = stiffness->at(O + K, R) = stiffness->at(O + L, P + M) = stiffness->at(O + M, P + L) = mass_coefs * m_val[N];
		}

		for(index_ts N = 0; N < current_stiffness.c_size; ++N) {
			const auto O = s_row[N], P = s_col[N], Q = O + K, R = O + L, S = P + K, T = P + L;
			stiffness->at(Q, T) = stiffness->at(R, S) = stiffness_coef * s_val[N];
			stiffness->at(O + J, S) = stiffness->at(Q, P + J) = stiffness->at(R, P + M) = stiffness->at(O + M, T) = stiffness_coefs * s_val[N];
		}

		I = current_pos;
		J = current_pos += n_block;
		K = current_pos += n_block;
		L = current_pos += n_block;
		M = current_pos += n_block;
		order -= 2.;
	}

	if(order < .5) assemble_by_mode_zero(current_pos = J, mass_coef, stiffness_coef);
	else {
		// eq. 53 
		for(index_tm N = 0; N < current_mass.c_size; ++N) {
			const auto O = m_row[N], P = m_col[N], Q = O + J, R = P + J;
			stiffness->at(O + L, P + L) = -(stiffness->at(O + I, R) = stiffness->at(Q, P + I) = mass_coef * m_val[N]);
			stiffness->at(Q, P + K) = stiffness->at(O + K, R) = mass_coefs * m_val[N];
		}

		for(index_ts N = 0; N < current_stiffness.c_size; ++N) {
			const auto O = s_row[N], P = s_col[N], Q = O + K, R = P + K, S = O + L, T = P + L;
			stiffness->at(S, T) = -(stiffness->at(Q, T) = stiffness->at(S, R) = stiffness_coef * s_val[N]);
			stiffness->at(O + J, R) = stiffness->at(Q, P + J) = stiffness_coefs * s_val[N];
		}
	}
}

void LeeNewmarkFull::assemble_by_mode_two(uword& current_pos, double mass_coef, double stiffness_coef, const double left_order, const double right_order) const {
	const auto s_order = left_order + right_order + 1.;
	const auto r = (2. * left_order + 1.) / (2. * right_order + 1.);
	const auto a = .5 * (1. + r) * pow(r, (-1. - left_order) / s_order);
	const auto b = a * pow(r, 1. / s_order);

	mass_coef *= a;
	stiffness_coef *= b;

	const auto& m_row = current_mass.row_mem();
	const auto& m_col = current_mass.col_mem();
	const auto& m_val = current_mass.val_mem();
	const auto& s_row = current_stiffness.row_mem();
	const auto& s_col = current_stiffness.col_mem();
	const auto& s_val = current_stiffness.val_mem();

	auto I = current_pos, K = current_pos;
	auto J = current_pos + n_block * static_cast<uword>(std::max(1., right_order));

	for(index_tm L = 0; L < current_mass.c_size; ++L) {
		const auto M = m_row[L], N = m_col[L], O = M + I, P = N + I;
		stiffness->at(O, N) = stiffness->at(O, N + J) = stiffness->at(M, P) = stiffness->at(M + J, P) = mass_coef * m_val[L];
	}

	auto formulate_block = [&](double order, double m_coef, double s_coef, const bool sign) {
		I = current_pos;
		J = current_pos += n_block;
		K = current_pos += n_block;

		while(order > 1.5) {
			for(index_tm L = 0; L < current_mass.c_size; ++L) {
				const auto M = m_row[L], N = m_col[L];
				stiffness->at(N + J, M + K) = stiffness->at(N + K, M + J) = m_coef * m_val[L];
			}

			for(index_ts L = 0; L < current_stiffness.c_size; ++L) {
				const auto M = s_row[L], N = s_col[L];
				stiffness->at(M + I, N + J) = stiffness->at(M + J, N + I) = s_coef * s_val[L];
			}

			I = current_pos;
			J = current_pos += n_block;
			K = current_pos += n_block;
			order -= 2.;
		}

		if(sign) {
			m_coef = -m_coef;
			s_coef = -s_coef;
		}

		if(order > .5) {
			// eq. 68
			for(index_tm L = 0; L < current_mass.c_size; ++L) {
				const auto M = m_row[L], N = m_col[L];
				stiffness->at(N + J, M + J) = m_coef * m_val[L];
			}
			for(index_ts L = 0; L < current_stiffness.c_size; ++L) {
				const auto M = s_row[L], N = s_col[L];
				stiffness->at(M + I, N + J) = stiffness->at(M + J, N + I) = s_coef * s_val[L];
			}
			current_pos = K;
		}
		else if(order > -.5) {
			// eq. 73
			for(index_ts L = 0; L < current_stiffness.c_size; ++L) {
				const auto M = s_row[L], N = s_col[L];
				stiffness->at(M + I, N + I) = -s_coef * s_val[L];
			}
			current_pos = J;
		}
		else {
			// eq. 71
			for(index_tm L = 0; L < current_mass.c_size; ++L) {
				const auto M = m_row[L], N = m_col[L];
				stiffness->at(N + I, M + I) = -m_coef * m_val[L];
			}
			current_pos = J;
		}
	};

	// central block
	formulate_block(right_order - 1., mass_coef, stiffness_coef, false);

	// right bottom corner
	formulate_block(left_order, mass_coef, stiffness_coef, true);
}

void LeeNewmarkFull::assemble_by_mode_three(uword& current_pos, double mass_coef, double stiffness_coef, const double gm) const {
	mass_coef *= 1. + gm;      // eq. 30
	stiffness_coef *= 1. + gm; // eq. 30

	const auto mass_coefs = .25 / gm * mass_coef;                  // eq. 87
	const auto stiffness_coefs = (1. + .25 / gm) * stiffness_coef; // eq. 87

	const auto I = current_pos;
	const auto J = current_pos += n_block;

	auto row = current_mass.row_mem();
	auto col = current_mass.col_mem();
	auto val = current_mass.val_mem();
	for(index_tm L = 0; L < current_mass.c_size; ++L) {
		const auto M = row[L], N = col[L], O = M + I, P = M + J, R = N + I, S = N + J;
		stiffness->at(O, N) = stiffness->at(M, R) = -(stiffness->at(M, N) = stiffness->at(O, R) = mass_coef * val[L]);
		stiffness->at(P, S) = mass_coefs * val[L];
	}
	row = current_stiffness.row_mem();
	col = current_stiffness.col_mem();
	val = current_stiffness.val_mem();
	for(index_ts L = 0; L < current_stiffness.c_size; ++L) {
		const auto M = row[L], N = col[L], O = M + I, P = M + J, R = N + I, S = N + J;
		stiffness->at(P, R) = stiffness->at(O, S) = -(stiffness->at(O, R) = stiffness_coef * val[L]);
		stiffness->at(P, S) = stiffness_coefs * val[L];
	}

	current_pos += n_block;
}

LeeNewmarkFull::LeeNewmarkFull(const unsigned T, std::vector<Mode>&& M, const double A, const double B)
	: LeeNewmarkBase(T, A, B)
	, damping_mode(std::forward<std::vector<Mode>>(M)) { for(auto& I : damping_mode) if(Type::T3 == I.t && (I.p.empty() || I.p.front() == 0. || I.p.front() < -1.)) I.t = Type::T0; }

int LeeNewmarkFull::initialize() {
	if(SUANPAN_SUCCESS != LeeNewmarkBase::initialize()) return SUANPAN_FAIL;

	if(StorageScheme::SPARSE != factory->get_storage_scheme()) {
		suanpan_error("please use command `set sparse_mat true` to enable sparse storage.\n");
		return SUANPAN_FAIL;
	}

	return SUANPAN_SUCCESS;
}

int LeeNewmarkFull::process_constraint() {
	const auto& D = get_domain().lock();

	// process constraint for the first time to obtain proper stiffness
	if(SUANPAN_SUCCESS != Integrator::process_constraint()) return SUANPAN_FAIL;

	auto& t_stiff = factory->get_stiffness()->triplet_mat;

	t_stiff.csc_condense();

	if(first_iteration) {
		// preallocate memory
		stiffness->triplet_mat.resize(get_amplifier() * t_stiff.c_size);
		stiffness->zeros();

		// ! deal with mass matrix first
		// the intact mass matrix will be the correct mass to be used
		// D->assemble_current_mass();

		auto& t_mass = factory->get_mass()->triplet_mat;

		access::rw(current_mass) = std::move(t_mass);

		// must initialise it since nothing will be checked in if left uninitialised
		t_mass = triplet_form<double, uword>(n_block, n_block, current_mass.c_size);

		// ! now deal with stiffness matrix
		auto& t_rabbit = access::rw(rabbit);

		// use local variable to temporarily store the original effective stiffness matrix
		t_rabbit = std::move(t_stiff);

		// preallocate to formulate current stiffness matrix
		t_stiff = triplet_form<double, uword>(n_block, n_block, t_rabbit.c_size);

		D->assemble_current_stiffness();
		if(SUANPAN_SUCCESS != Integrator::process_constraint()) return SUANPAN_FAIL;
		t_stiff.csc_condense();

		access::rw(current_stiffness) = std::move(t_stiff);

		// now current mass and stiffness are formulated
		// assemble unrolled damping matrix and the corresponding damping force
		update_stiffness();
		update_residual();

		// move original effective stiffness matrix back
		// need to add it to global stiffness matrix later
		t_stiff = std::move(t_rabbit);

		const auto& row = stiffness->triplet_mat.row_mem();
		const auto& col = stiffness->triplet_mat.col_mem();
		const auto& val = stiffness->triplet_mat.val_mem();

		t_rabbit = triplet_form<double, uword>(n_block, n_block, t_stiff.c_size);
		for(decltype(stiffness->triplet_mat)::index_type I = 0; I < stiffness->triplet_mat.c_size; ++I) {
			// quit if current column is beyond the original size of matrix
			if(col[I] >= n_block) break;
			// check in left top block of unrolled damping matrix to be used in subsequent iterations
			if(row[I] < n_block) t_rabbit.at(row[I], col[I]) = val[I];
		}

		stiffness->triplet_mat += t_stiff;

		first_iteration = false;
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

		// check in original nonzero entries in unrolled damping matrix
		stiffness->triplet_mat += rabbit;

		update_residual();

		stiffness->triplet_mat += t_stiff;
	}

	return SUANPAN_SUCCESS;
}

void LeeNewmarkFull::print() { suanpan_info("A Newmark solver using Lee's damping model with adjustable bandwidth. DOI: 10.1016/j.compstruc.2020.106423\n"); }
