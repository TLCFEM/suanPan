/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "LeeNewmarkIterative.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

void LeeNewmarkIterative::init_worker(const unsigned n_dim, const unsigned n_multiplier) {
    auto& W = get_domain()->get_factory();

    const auto [l,u] = W->get_bandwidth();

    worker = make_unique<SparseMatSuperLU<double>>(n_dim, n_dim, n_multiplier * (l + u + 1) * W->get_size());
}

void LeeNewmarkIterative::assemble(const shared_ptr<MetaMat<double>>& in_mat, const uword row_shift, const uword col_shift, const double scalar) const {
    auto& W = get_domain()->get_factory();

    const auto n_block = W->get_size();

    if(W->is_sparse()) worker->triplet_mat.assemble(in_mat->triplet_mat, row_shift, col_shift, scalar);
    else {
        const auto [low, up] = W->get_bandwidth();
        for(unsigned L = 0; L < n_block; ++L) {
            const auto N = L + col_shift;
            for(unsigned K = std::max(L, up) - up; K < std::min(n_block, L + low + 1); ++K) {
                const auto M = K + row_shift;
                if(const sp_d auto t_val = in_mat->operator()(K, L); !suanpan::approx_equal(0., t_val)) worker->at(M, N) = scalar * t_val;
            }
        }
    }
}

void LeeNewmarkIterative::assemble_mass(const uword row_shift, const uword col_shift, const double scalar) const { assemble(current_mass, row_shift, col_shift, scalar); }

void LeeNewmarkIterative::assemble_stiffness(const uword row_shift, const uword col_shift, const double scalar) const { assemble(current_stiffness, row_shift, col_shift, scalar); }

void LeeNewmarkIterative::assemble_mass(const std::vector<sword>& row_shift, const std::vector<sword>& col_shift, const std::vector<double>& scalar) const {
    suanpan_assert([&] { if(scalar.size() != row_shift.size() || scalar.size() != col_shift.size()) throw invalid_argument("size mismatch detected"); });

    for(decltype(scalar.size()) I = 0; I < scalar.size(); ++I) if(row_shift[I] >= 0 && col_shift[I] >= 0) assemble_mass(row_shift[I], col_shift[I], scalar[I]);
}

void LeeNewmarkIterative::assemble_stiffness(const std::vector<sword>& row_shift, const std::vector<sword>& col_shift, const std::vector<double>& scalar) const {
    suanpan_assert([&] { if(scalar.size() != row_shift.size() || scalar.size() != col_shift.size()) throw invalid_argument("size mismatch detected"); });

    for(decltype(scalar.size()) I = 0; I < scalar.size(); ++I) if(row_shift[I] >= 0 && col_shift[I] >= 0) assemble_stiffness(row_shift[I], col_shift[I], scalar[I]);
}

vec LeeNewmarkIterative::update_by_mode_zero(const double mass_coef, const double stiffness_coef) const {
    auto& W = get_domain()->get_factory();

    const auto kernel = current_mass->make_copy();
    kernel += stiffness_coef / mass_coef * current_stiffness;
    const auto damping_force = current_mass * W->get_trial_velocity();
    vec tmp;
    kernel->solve(tmp, damping_force);
    return mass_coef * (damping_force - current_mass * tmp);
}

vec LeeNewmarkIterative::update_by_mode_one(const double mass_coef, const double stiffness_coef, int order) {
    auto& W = get_domain()->get_factory();

    const int64_t n_block = W->get_size();
    const auto n_total = (2 * order + 1) * n_block;

    init_worker(static_cast<unsigned>(n_total), 3llu + 6llu * order);

    const auto mass_coefs = .5 * mass_coef;           // eq. 10
    const auto stiffness_coefs = .5 * stiffness_coef; // eq. 10

    auto current_pos = 0ll;

    auto I = -n_block;
    auto J = current_pos;
    auto K = current_pos += n_block;
    auto L = current_pos += n_block;
    auto M = current_pos += n_block;

    while(order > 1) {
        // eq. 61
        assemble_mass({I, J, J, K, L, M}, {J, I, K, J, M, L}, {mass_coef, mass_coef, mass_coefs, mass_coefs, mass_coefs, mass_coefs});
        assemble_stiffness({K, L, J, K, L, M}, {L, K, K, J, M, L}, {stiffness_coef, stiffness_coef, stiffness_coefs, stiffness_coefs, stiffness_coefs, stiffness_coefs});

        I = current_pos;
        J = current_pos += n_block;
        K = current_pos += n_block;
        L = current_pos += n_block;
        M = current_pos += n_block;
        order -= 2;
    }

    // eq. 3
    if(order < 1) {
        assemble_mass({I, J, I, J}, {I, J, J, I}, {mass_coef, mass_coef, mass_coef, mass_coef});
        assemble_stiffness(J, J, stiffness_coef);
    }
    else {
        // eq. 53
        assemble_mass({I, J, J, K, L}, {J, I, K, J, L}, {mass_coef, mass_coef, mass_coefs, mass_coefs, -mass_coef});
        assemble_stiffness({J, K, K, L, L}, {K, J, L, K, L}, {stiffness_coefs, stiffness_coefs, stiffness_coef, stiffness_coef, -stiffness_coef});
    }

    vec damping_force(n_total, fill::zeros);
    damping_force.head(n_block) = current_mass * W->get_trial_velocity();
    vec tmp_a;
    worker->solve(tmp_a, damping_force);
    const vec tmp_b = -mass_coef * mass_coef * tmp_a.head(n_block);
    return current_mass * tmp_b;
}

vec LeeNewmarkIterative::update_by_mode_two(double, double, int, int) const { return {}; }

vec LeeNewmarkIterative::update_by_mode_three(double mass_coef, double stiffness_coef, const double gm) {
    auto& W = get_domain()->get_factory();

    const auto n_block = W->get_size();
    const auto n_total = 2 * n_block;

    init_worker(n_total, 6);

    mass_coef *= 1. + gm;      // eq. 30
    stiffness_coef *= 1. + gm; // eq. 30

    constexpr auto I = 0;
    const auto J = n_block;

    assemble_mass({J, I}, {J, I}, {.25 / gm * mass_coef, mass_coef});
    assemble_stiffness({J, I, I, J}, {J, I, J, I}, {(1. + .25 / gm) * stiffness_coef, stiffness_coef, -stiffness_coef, -stiffness_coef});

    vec damping_force(n_total, fill::zeros);
    damping_force.head(n_block) = current_mass * W->get_trial_velocity() * mass_coef;
    vec tmp_a;
    worker->solve(tmp_a, damping_force);
    const vec tmp_b = mass_coef * tmp_a.head(n_block);
    return damping_force.head(n_block) - current_mass * tmp_b;
}

vec LeeNewmarkIterative::update_by_mode_four(double, double, int, int, int, int, double) const { return {}; }

void LeeNewmarkIterative::update_damping_force() {
    auto& W = get_domain()->get_factory();

    vec summation(W->get_size(), fill::zeros);

    const auto i = [](const double x) { return static_cast<int>(x); };

    for(const auto& [t, p, zeta, omega] : damping_mode) {
        const auto mass_coef = 4. * zeta * omega, stiffness_coef = 4. * zeta / omega;
        switch(t) {
        case Type::T0:
            summation += update_by_mode_zero(mass_coef, stiffness_coef);
            break;
        case Type::T1:
            summation += update_by_mode_one(mass_coef, stiffness_coef, i(p.front()));
            break;
        case Type::T2:
            summation += update_by_mode_two(mass_coef, stiffness_coef, i(p(0)), i(p(1)));
            break;
        case Type::T3:
            summation += update_by_mode_three(mass_coef, stiffness_coef, p.front());
            break;
        case Type::T4:
            summation += update_by_mode_four(mass_coef, stiffness_coef, i(p(0)), i(p(1)), i(p(2)), i(p(3)), p(4));
            break;
        }
    }

    W->update_trial_damping_force_by(summation);
    W->update_sushi_by(summation);
}

LeeNewmarkIterative::LeeNewmarkIterative(const unsigned T, std::vector<Mode>&& M, const double A, const double B)
    : Newmark(T, A, B)
    , damping_mode(std::move(M)) {
    for(auto& [t, p, zeta, omega] : damping_mode)
        switch(t) {
        case Type::T0:
            break;
        case Type::T1:
            if(suanpan::approx_equal(p(0), 0.)) t = Type::T0;
            break;
        case Type::T2:
            if(suanpan::approx_equal(p(0) + p(1), 0.)) t = Type::T0;
            break;
        case Type::T3:
            if(suanpan::approx_equal(p(0), 0.) || p(0) < -1.) t = Type::T0;
            break;
        case Type::T4:
            if(p(4) < -1.) t = Type::T0;
            else if(suanpan::approx_equal(p(4), 0.)) t = suanpan::approx_equal(p(0) + p(1), 0.) ? Type::T0 : Type::T2;
            else if(suanpan::approx_equal(p(0) + p(1) + p(2) + p(3), 0.)) {
                t = Type::T3;
                p = vec{p(4)};
            }
            break;
        }
}

int LeeNewmarkIterative::process_constraint() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    auto& t_mass = W->modify_mass();
    auto& t_stiffness = W->modify_stiffness();

    current_mass.swap(t_mass);
    current_stiffness.swap(t_stiffness);
    if(SUANPAN_SUCCESS != Newmark::process_constraint()) return SUANPAN_FAIL;
    current_mass.swap(t_mass);
    current_stiffness.swap(t_stiffness);

    update_damping_force();

    return Newmark::process_constraint();
}

int LeeNewmarkIterative::process_constraint_resistance() {
    update_damping_force();

    return Newmark::process_constraint_resistance();
}

void LeeNewmarkIterative::assemble_matrix() {
    const auto D = get_domain();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_trial_stiffness(); });
    auto fb = std::async([&] { D->assemble_trial_geometry(); });
    auto fc = std::async([&] { D->assemble_trial_damping(); });
    auto fd = std::async([&] { D->assemble_trial_nonviscous(); });
    auto fe = std::async([&] { D->assemble_trial_mass(); });

    fa.get();
    fb.get();
    fc.get();
    fd.get();
    fe.get();

    if(W->is_nlgeom()) W->get_stiffness() += W->get_geometry();

    current_mass = W->get_mass()->make_copy();
    current_stiffness = W->get_stiffness()->make_copy();

    W->get_stiffness() += C0 * W->get_mass();

    W->get_stiffness() += W->is_nonviscous() ? C1 * (W->get_damping() + W->get_nonviscous()) : C1 * W->get_damping();
}
