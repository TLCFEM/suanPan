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
    const auto [l,u] = factory->get_bandwidth();

    if(SolverType::MUMPS == factory->get_solver_type()) worker = make_unique<SparseMatMUMPS<double>>(n_dim, n_dim, n_multiplier * (l + u + 1) * n_block);
    else if(SolverType::LIS == factory->get_solver_type()) worker = make_unique<SparseMatLis<double>>(n_dim, n_dim, n_multiplier * (l + u + 1) * n_block);
#ifdef SUANPAN_MKL
    else if(SolverType::PARDISO == factory->get_solver_type()) worker = make_unique<SparseMatPARDISO<double>>(n_dim, n_dim, n_multiplier * (l + u + 1) * n_block);
    else if(SolverType::FGMRES == factory->get_solver_type()) worker = make_unique<SparseMatFGMRES<double>>(n_dim, n_dim, n_multiplier * (l + u + 1) * n_block);
#endif
#ifdef SUANPAN_CUDA
    else if(SolverType::CUDA == factory->get_solver_type()) worker = make_unique<SparseMatCUDA<double>>(n_dim, n_dim, n_multiplier * (l + u + 1) * n_block);
#endif
    else worker = make_unique<SparseMatSuperLU<double>>(n_dim, n_dim, n_multiplier * (l + u + 1) * n_block);
}

void LeeNewmarkIterative::assemble(const shared_ptr<MetaMat<double>>& in_mat, const uword row_shift, const uword col_shift, const double scalar) const {
    if(factory->is_sparse()) worker->triplet_mat.assemble(in_mat->triplet_mat, row_shift, col_shift, scalar);
    else {
        const auto [low, up] = factory->get_bandwidth();
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

void LeeNewmarkIterative::formulate_block(sword& current_pos, const double m_coef, const double s_coef, int order) const {
    auto I = current_pos;
    auto J = current_pos += n_block;
    auto K = current_pos += n_block;

    while(order > 1) {
        assemble_mass({J, K}, {K, J}, {m_coef, m_coef});
        assemble_stiffness({I, J}, {J, I}, {s_coef, s_coef});

        I = current_pos;
        J = current_pos += n_block;
        K = current_pos += n_block;
        order -= 2;
    }

    if(order > 0) {
        // eq. 68
        assemble_mass(J, J, -m_coef);
        assemble_stiffness({I, J}, {J, I}, {s_coef, s_coef});

        current_pos = K;

        return;
    }

    if(order > -1) {
        // eq. 73
        assemble_stiffness(I, I, s_coef);

        current_pos = J;

        return;
    }

    // eq. 71
    assemble_mass(I, I, m_coef);

    current_pos = J;
}

void LeeNewmarkIterative::formulate_block(sword& current_pos, const std::vector<double>& m_coef, const std::vector<double>& s_coef, const std::vector<int>& order) const {
    suanpan_assert([&] { if(order.size() != m_coef.size() || order.size() != s_coef.size()) throw invalid_argument("size mismatch detected"); });

    for(size_t I = 0; I < order.size(); ++I) formulate_block(current_pos, m_coef[I], s_coef[I], order[I]);
}

vec LeeNewmarkIterative::update_by_mode_zero(const double mass_coef, const double stiffness_coef) const {
    const auto kernel = current_mass->make_copy();
    kernel += stiffness_coef / mass_coef * current_stiffness;
    const auto damping_force = current_mass * factory->get_trial_velocity();
    return mass_coef * (damping_force - current_mass * kernel->solve(damping_force));
}

vec LeeNewmarkIterative::update_by_mode_one(const double mass_coef, const double stiffness_coef, int order) {
    const auto n_total = (2 * order + 1) * n_block;

    init_worker(n_total, 3llu + 6llu * order);

    const auto mass_coefs = .5 * mass_coef;           // eq. 10
    const auto stiffness_coefs = .5 * stiffness_coef; // eq. 10

    sword current_pos{0};

    auto I = -static_cast<sword>(n_block);
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
    damping_force.head(n_block) = current_mass * factory->get_trial_velocity();
    const vec tmp_b = -mass_coef * mass_coef * worker->solve(damping_force).head_rows(n_block);
    return current_mass * tmp_b;
}

vec LeeNewmarkIterative::update_by_mode_two(double mass_coef, double stiffness_coef, const int npr, const int npl) {
    const auto n_total = (npr + npl + 1) * n_block;

    init_worker(n_total, 2u + 5u * static_cast<unsigned>(.5 * (npr + npl - 1.)));

    const auto nps = npr + npl + 1.;
    const auto r = (2. * npl + 1.) / (2. * npr + 1.);
    const auto a = .5 * (1. + r) * pow(r, (-1. - npl) / nps);
    const auto b = a * pow(r, 1. / nps);

    mass_coef *= a;      // eq. 18
    stiffness_coef *= b; // eq. 18

    sword current_pos{0};

    vec damping_force(n_total, fill::zeros), final_force;

    // eq. 100
    if(0 == npr) {
        assemble_mass(0, 0, mass_coef);

        formulate_block(current_pos, mass_coef, stiffness_coef, npl);

        damping_force.head(n_block) = current_mass * factory->get_trial_velocity() * mass_coef;
        const vec tmp_b = mass_coef * worker->solve(damping_force).head_rows(n_block);
        final_force = damping_force.head(n_block) - current_mass * tmp_b;
    }
    else {
        const auto J = n_block * npr;

        // eq. 81
        assemble_mass({0, J}, {J, 0}, {mass_coef, mass_coef});

        // central block, right bottom corner
        formulate_block(current_pos, {-mass_coef, mass_coef}, {-stiffness_coef, stiffness_coef}, {npr - 1, npl});

        damping_force.head(n_block) = current_mass * factory->get_trial_velocity();
        const vec tmp_b = -mass_coef * mass_coef * worker->solve(damping_force).head_rows(n_block);
        final_force = current_mass * tmp_b;
    }

    return final_force;
}

vec LeeNewmarkIterative::update_by_mode_three(double mass_coef, double stiffness_coef, const double gm) {
    const auto n_total = 2 * n_block;

    init_worker(n_total, 6);

    mass_coef *= 1. + gm;      // eq. 30
    stiffness_coef *= 1. + gm; // eq. 30

    constexpr auto I = 0;
    const auto J = n_block;

    assemble_mass({J, I}, {J, I}, {.25 / gm * mass_coef, mass_coef});
    assemble_stiffness({J, I, I, J}, {J, I, J, I}, {(1. + .25 / gm) * stiffness_coef, stiffness_coef, -stiffness_coef, -stiffness_coef});

    vec damping_force(n_total, fill::zeros);
    damping_force.head(n_block) = current_mass * factory->get_trial_velocity() * mass_coef;
    const vec tmp_b = mass_coef * worker->solve(damping_force).head_rows(n_block);
    return damping_force.head(n_block) - current_mass * tmp_b;
}

vec LeeNewmarkIterative::update_by_mode_four(const double mass_coef, const double stiffness_coef, const int npr, const int npl, const int npk, const int npm, const double gm) {
    const auto n_total = (npr + npl + npk + npm + 2) * n_block;

    init_worker(n_total, 2u * (npr + npl + npk + npm) + 8u);

    const auto rs = (2. * npl + 1.) / (2. * npr + 1.);
    const auto rp = (2. * npm + 1.) / (2. * npk + 1.);
    const auto nps = npr + npl + 1.;
    const auto npt = npk + npm + 1.;
    const auto as = 2. * (1. + rs) * pow(rs, (-1. - npl) / nps);
    const auto ap = 2. * (1. + rp) * pow(rp, (-1. - npm) / npt);
    const auto bs = as * pow(rs, 1. / nps);
    const auto bp = ap * pow(rp, 1. / npt);
    const auto bgm = .25 * ap * bp * gm;

    const auto m_coef_s = .25 * (1. + gm) * as * mass_coef;
    const auto m_coef_p = .25 * (1. + gm) * ap * mass_coef;
    const auto s_coef_s = .25 * (1. + gm) * bs * stiffness_coef;
    const auto s_coef_p = .25 * (1. + gm) * bp * stiffness_coef;

    sword current_pos{0};

    vec damping_force(n_total, fill::zeros), final_force;
    damping_force.head(n_block) = current_mass * factory->get_trial_velocity() * m_coef_s;

    const auto solve_a = [&](const uword middle) {
        damping_force.subvec(middle, middle + n_block) = damping_force.head(n_block);
        const vec tmp_a = worker->solve(damping_force);
        const vec tmp_b = m_coef_s * (tmp_a.head(n_block) + tmp_a.subvec(middle, middle + n_block));
        final_force = damping_force.head(n_block) - current_mass * tmp_b;
    };

    const auto solve_b = [&] {
        const vec tmp_b = -m_coef_s * worker->solve(damping_force).head_rows(n_block);
        final_force = current_mass * tmp_b;
    };

    if(0 == npr && 0 == npm) {
        // eq. 100
        const auto J = static_cast<sword>(n_block) * npl + n_block;

        assemble_mass({0, 0, J, J}, {0, J, 0, J}, {m_coef_s, m_coef_s, m_coef_s, m_coef_s + m_coef_p / bgm});

        formulate_block(current_pos, {m_coef_s, m_coef_p / bgm}, {s_coef_s, s_coef_p / bgm}, {npl, npk});

        solve_a(J);
    }
    else if(0 == npr) {
        // eq. 98
        const auto J = static_cast<sword>(n_block) * npl + n_block;
        const auto K = J + static_cast<sword>(n_block) * npk + n_block;

        assemble_mass({0, 0, J, J, J, K}, {0, J, 0, J, K, J}, {m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_p, m_coef_p});

        formulate_block(current_pos, {m_coef_s, m_coef_p / bgm, -bgm * m_coef_p}, {s_coef_s, s_coef_p / bgm, -bgm * s_coef_p}, {npl, npk, npm - 1});

        solve_a(J);
    }
    else if(0 == npm) {
        // eq. 97
        const auto J = static_cast<sword>(n_block) * npr;
        const auto K = J + static_cast<sword>(n_block) * npl + n_block;

        assemble_mass({0, 0, J, K, K}, {J, K, 0, 0, K}, {m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_p / bgm});

        formulate_block(current_pos, {-m_coef_s, m_coef_s, m_coef_p / bgm}, {-s_coef_s, s_coef_s, s_coef_p / bgm}, {npr - 1, npl, npk});

        solve_b();
    }
    else {
        // eq. 84
        const auto J = static_cast<sword>(n_block) * npr;
        const auto K = J + static_cast<sword>(n_block) * npl + n_block;
        const auto L = K + static_cast<sword>(n_block) * npk + n_block;

        assemble_mass({0, 0, J, K, K, L}, {J, K, 0, 0, L, K}, {m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_p, m_coef_p});

        formulate_block(current_pos, {-m_coef_s, m_coef_s, m_coef_p / bgm, -bgm * m_coef_p}, {-s_coef_s, s_coef_s, s_coef_p / bgm, -bgm * s_coef_p}, {npr - 1, npl, npk, npm - 1});

        solve_b();
    }

    return final_force;
}

void LeeNewmarkIterative::update_damping_force() {
    vec summation(n_block, fill::zeros);

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

    factory->update_trial_damping_force_by(summation);
    factory->update_sushi_by(summation);
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

int LeeNewmarkIterative::initialize() {
    if(Newmark::initialize() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    factory = get_domain()->get_factory();

    access::rw(n_block) = factory->get_size();

    return SUANPAN_SUCCESS;
}

int LeeNewmarkIterative::process_constraint() {
    auto& t_mass = factory->modify_mass();
    auto& t_stiffness = factory->modify_stiffness();

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

void LeeNewmarkIterative::print() {
    suanpan_info("A Newmark solver using Lee's damping model with iterative solving strategy.\n");
}
