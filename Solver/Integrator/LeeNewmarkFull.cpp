/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "LeeNewmarkFull.h"
#include <Domain/DomainBase.h>
#include <Domain/FactoryHelper.hpp>

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
        else if(Type::T4 == t) n_size += 2llu * static_cast<uword>(p(0) + p(1) + p(2) + p(3) + 4.);

    return n_size;
}

/**
 * \brief compute the exact size of final global effective stiffness
 * \return the exact size of final global effective stiffness
 */
uword LeeNewmarkFull::get_total_size() const {
    auto n_size = 1llu;

    for(const auto& [t, p, zeta, omega] : damping_mode)
        if(Type::T0 == t) n_size += 1llu;
        else if(Type::T1 == t) n_size += 2llu * static_cast<uword>(p.front()) + 1llu;
        else if(Type::T2 == t) n_size += static_cast<uword>(p(0) + p(1)) + 1llu;
        else if(Type::T3 == t) n_size += 2llu;
        else if(Type::T4 == t) n_size += static_cast<uword>(p(0) + p(1) + p(2) + p(3)) + 2llu;

    return n_size * n_block;
}

void LeeNewmarkFull::update_stiffness() const {
    if(build_graph) {
        const auto t_size = get_total_size() / n_block;
        access::rw(stiffness_graph).set_size(t_size, t_size);
        access::rw(mass_graph).set_size(t_size, t_size);
    }

    auto IDX = n_block;

    // ! make sure global stiffness only holds unrolled damping matrix when exit
    stiffness->zeros();

    for(const auto& [t, p, zeta, omega] : damping_mode)
        if(const auto mass_coef = 4. * zeta * omega * C1, stiffness_coef = 4. * zeta / omega * C1; Type::T0 == t) assemble_by_mode_zero(IDX, mass_coef, stiffness_coef);
        else if(Type::T1 == t) assemble_by_mode_one(IDX, mass_coef, stiffness_coef, static_cast<int>(p.front()));
        else if(Type::T2 == t) assemble_by_mode_two(IDX, mass_coef, stiffness_coef, static_cast<int>(p(0)), static_cast<int>(p(1)));
        else if(Type::T3 == t) assemble_by_mode_three(IDX, mass_coef, stiffness_coef, p.front());
        else if(Type::T4 == t) assemble_by_mode_four(IDX, mass_coef, stiffness_coef, static_cast<int>(p(0)), static_cast<int>(p(1)), static_cast<int>(p(2)), static_cast<int>(p(3)), p(4));

    stiffness->csc_condense();

    if(build_graph) access::rw(build_graph) = false;
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
    auto fa = std::async([&] { get_trial_damping_force(factory) -= residual.head(n_block); });
    auto fb = std::async([&] { get_incre_damping_force(factory) -= residual.head(n_block); });
    // ! update left-hand side
    auto fc = std::async([&] { get_sushi(factory) -= residual.head(n_block); });

    fa.get();
    fb.get();
    fc.get();
}

void LeeNewmarkFull::assemble_mass(const uword row_shift, const uword col_shift, const double scalar) const { assemble(access::rw(mass_graph), current_mass, row_shift, col_shift, scalar); }

void LeeNewmarkFull::assemble_stiffness(const uword row_shift, const uword col_shift, const double scalar) const { assemble(access::rw(stiffness_graph), current_stiffness, row_shift, col_shift, scalar); }

void LeeNewmarkFull::assemble_mass(const std::vector<uword>& row_shift, const std::vector<uword>& col_shift, const std::vector<double>& scalar) const {
    suanpan_debug([&] { if(scalar.size() != row_shift.size() || scalar.size() != col_shift.size()) throw invalid_argument("size mismatch detected"); });

    for(decltype(scalar.size()) I = 0; I < scalar.size(); ++I) assemble_mass(row_shift[I], col_shift[I], scalar[I]);
}

void LeeNewmarkFull::assemble_stiffness(const std::vector<uword>& row_shift, const std::vector<uword>& col_shift, const std::vector<double>& scalar) const {
    suanpan_debug([&] { if(scalar.size() != row_shift.size() || scalar.size() != col_shift.size()) throw invalid_argument("size mismatch detected"); });

    for(decltype(scalar.size()) I = 0; I < scalar.size(); ++I) assemble_stiffness(row_shift[I], col_shift[I], scalar[I]);
}

void LeeNewmarkFull::formulate_block(uword& current_pos, const double m_coef, const double s_coef, int order) const {
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

void LeeNewmarkFull::formulate_block(uword& current_pos, const std::vector<double>& m_coef, const std::vector<double>& s_coef, const std::vector<int>& order) const {
    suanpan_debug([&] { if(order.size() != m_coef.size() || order.size() != s_coef.size()) throw invalid_argument("size mismatch detected"); });

    for(size_t I = 0; I < order.size(); ++I) formulate_block(current_pos, m_coef[I], s_coef[I], order[I]);
}

void LeeNewmarkFull::assemble_by_mode_zero(uword& current_pos, const double mass_coef, const double stiffness_coef) const {
    const auto I = current_pos;

    assemble_mass({0, 0, I, I}, {0, I, 0, I}, {mass_coef, mass_coef, mass_coef, mass_coef});
    assemble_stiffness(I, I, stiffness_coef);

    current_pos += n_block;
}

void LeeNewmarkFull::assemble_by_mode_one(uword& current_pos, const double mass_coef, const double stiffness_coef, int order) const {
    const auto mass_coefs = .5 * mass_coef;           // eq. 10
    const auto stiffness_coefs = .5 * stiffness_coef; // eq. 10

    auto I = 0llu;
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

        current_pos = K;

        return;
    }

    // eq. 53
    assemble_mass({I, J, J, K, L}, {J, I, K, J, L}, {mass_coef, mass_coef, mass_coefs, mass_coefs, -mass_coef});
    assemble_stiffness({J, K, K, L, L}, {K, J, L, K, L}, {stiffness_coefs, stiffness_coefs, stiffness_coef, stiffness_coef, -stiffness_coef});
}

void LeeNewmarkFull::assemble_by_mode_two(uword& current_pos, double mass_coef, double stiffness_coef, const int npr, const int npl) const {
    const auto nps = static_cast<double>(npr) + static_cast<double>(npl) + 1.;
    const auto r = (2. * npl + 1.) / (2. * npr + 1.);
    const auto a = .5 * (1. + r) * pow(r, (-1. - npl) / nps);
    const auto b = a * pow(r, 1. / nps);

    mass_coef *= a;      // eq. 18
    stiffness_coef *= b; // eq. 18

    const auto I = current_pos;

    // eq. 100
    if(0 == npr) {
        assemble_mass({0, 0, I, I}, {0, I, 0, I}, {mass_coef, mass_coef, mass_coef, mass_coef});

        formulate_block(current_pos, mass_coef, stiffness_coef, npl);

        return;
    }

    const auto J = current_pos + n_block * npr;

    // eq. 81
    assemble_mass({0, I, I, J}, {I, 0, J, I}, {mass_coef, mass_coef, mass_coef, mass_coef});

    // central block, right bottom corner
    formulate_block(current_pos, {-mass_coef, mass_coef}, {-stiffness_coef, stiffness_coef}, {npr - 1, npl});
}

void LeeNewmarkFull::assemble_by_mode_three(uword& current_pos, double mass_coef, double stiffness_coef, const double gm) const {
    mass_coef *= 1. + gm;      // eq. 30
    stiffness_coef *= 1. + gm; // eq. 30

    const auto I = current_pos;
    const auto J = current_pos += n_block;

    assemble_mass({J, 0, I, 0, I}, {J, 0, I, I, 0}, {.25 / gm * mass_coef, mass_coef, mass_coef, -mass_coef, -mass_coef});
    assemble_stiffness({J, I, I, J}, {J, I, J, I}, {(1. + .25 / gm) * stiffness_coef, stiffness_coef, -stiffness_coef, -stiffness_coef});

    current_pos += n_block;
}

void LeeNewmarkFull::assemble_by_mode_four(uword& current_pos, const double mass_coef, const double stiffness_coef, const int npr, const int npl, const int npk, const int npm, const double gm) const {
    const auto rs = (2. * npl + 1.) / (2. * npr + 1.);
    const auto rp = (2. * npm + 1.) / (2. * npk + 1.);
    const auto nps = static_cast<double>(npr) + static_cast<double>(npl) + 1.;
    const auto npt = static_cast<double>(npk) + static_cast<double>(npm) + 1.;
    const auto as = 2. * (1. + rs) * pow(rs, (-1. - npl) / nps);
    const auto ap = 2. * (1. + rp) * pow(rp, (-1. - npm) / npt);
    const auto bs = as * pow(rs, 1. / nps);
    const auto bp = ap * pow(rp, 1. / npt);
    const auto bgm = .25 * ap * bp * gm;

    const auto m_coef_s = .25 * (1. + gm) * as * mass_coef;
    const auto m_coef_p = .25 * (1. + gm) * ap * mass_coef;
    const auto s_coef_s = .25 * (1. + gm) * bs * stiffness_coef;
    const auto s_coef_p = .25 * (1. + gm) * bp * stiffness_coef;

    // eq. 100
    if(0 == npr && 0 == npm) {
        const auto I = current_pos;
        const auto J = I + n_block * npl + n_block;

        assemble_mass({0, 0, 0, I, I, I, J, J, J}, {0, I, J, 0, I, J, 0, I, J}, {m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s + m_coef_p / bgm});

        formulate_block(current_pos, {m_coef_s, m_coef_p / bgm}, {s_coef_s, s_coef_p / bgm}, {npl, npk});

        return;
    }

    // eq. 98
    if(0 == npr) {
        const auto I = current_pos;
        const auto J = I + n_block * npl + n_block;
        const auto K = J + n_block * npk + n_block;

        assemble_mass({0, 0, 0, I, I, I, J, J, J, J, K}, {0, I, J, 0, I, J, 0, I, J, K, J}, {m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_p, m_coef_p});

        formulate_block(current_pos, {m_coef_s, m_coef_p / bgm, -bgm * m_coef_p}, {s_coef_s, s_coef_p / bgm, -bgm * s_coef_p}, {npl, npk, npm - 1});

        return;
    }

    // eq. 97
    if(0 == npm) {
        const auto I = current_pos;
        const auto J = I + n_block * npr;
        const auto K = J + n_block * npl + n_block;

        assemble_mass({0, I, I, I, J, K, K}, {I, 0, J, K, I, I, K}, {m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_p / bgm});

        formulate_block(current_pos, {-m_coef_s, m_coef_s, m_coef_p / bgm}, {-s_coef_s, s_coef_s, s_coef_p / bgm}, {npr - 1, npl, npk});

        return;
    }

    // eq. 84
    const auto I = current_pos;
    const auto J = I + n_block * npr;
    const auto K = J + n_block * npl + n_block;
    const auto L = K + n_block * npk + n_block;

    assemble_mass({0, I, I, I, J, K, K, L}, {I, 0, J, K, I, I, L, K}, {m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_s, m_coef_p, m_coef_p});

    formulate_block(current_pos, {-m_coef_s, m_coef_s, m_coef_p / bgm, -bgm * m_coef_p}, {-s_coef_s, s_coef_s, s_coef_p / bgm, -bgm * s_coef_p}, {npr - 1, npl, npk, npm - 1});
}

LeeNewmarkFull::LeeNewmarkFull(const unsigned T, std::vector<Mode>&& M, const double A, const double B, const StiffnessType ST)
    : LeeNewmarkBase(T, A, B, ST)
    , damping_mode(std::forward<std::vector<Mode>>(M)) {
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

int LeeNewmarkFull::initialize() {
    if(SUANPAN_SUCCESS != LeeNewmarkBase::initialize()) return SUANPAN_FAIL;

    if(if_iterative && PreconditionerType::ILU != factory->get_solver_setting().preconditioner_type) {
        suanpan_error("iterative solver with preconditioner other than ILU is not supported, please consider LeeNewmark.\n");
        return SUANPAN_FAIL;
    }

    if(factory->is_sparse()) return SUANPAN_SUCCESS;

    suanpan_error("please use command `set sparse_mat true` to enable sparse storage.\n");
    return SUANPAN_FAIL;
}

int LeeNewmarkFull::process_constraint() {
    const auto& D = get_domain();

    // process constraint for the first time to obtain proper stiffness
    if(SUANPAN_SUCCESS != LeeNewmarkBase::process_constraint()) return SUANPAN_FAIL;

    // this stiffness contains geometry, mass and damping from Newmark::assemble_matrix()
    auto& t_stiff = factory->get_stiffness()->triplet_mat;

    t_stiff.csc_condense();

    const auto num_entry = 2 * t_stiff.n_elem;

    auto& t_triplet = stiffness->triplet_mat;

    // global matrix needs to be assembled as long as it is the first iteration
    // or trial stiffness matrix is used
    if(first_iteration || StiffnessType::TRIAL == stiffness_type) {
        // preallocate memory
        t_triplet.init(get_amplifier() * num_entry);
        // stiffness->zeros();

        // ! deal with mass matrix first
        // the intact mass matrix will be the correct mass to be used
        // assemble current mass if mass is changing
        // D->assemble_current_mass();

        // const auto& t_mass = factory->get_mass()->triplet_mat;

        // otherwise, directly make a copy
        auto f_mass = std::async([&] {
            auto& t_mass = factory->get_mass()->triplet_mat;
            t_mass.csc_condense();
            access::rw(current_mass) = t_mass;
        });

        // ! probably faster but unsafe
        // steal the already assembled mass matrix
        // access::rw(current_mass) = std::move(t_mass);
        // must initialise it since nothing will be checked in if left uninitialized
        // t_mass = triplet_form<double, uword>(n_block, n_block, current_mass.n_elem);

        // handle geometry matrix
        auto f_geometry = std::async([&] {
            using mat_t = triplet_form<double, uword>;

            if(!factory->is_nlgeom()) return mat_t();

            auto& t_geometry = factory->get_geometry()->triplet_mat;

            mat_t t_fox;

            // backup whatever
            std::swap(t_fox, t_geometry);

            // now t_geometry is empty
            t_geometry = triplet_form<double, uword>(n_block, n_block, num_entry);

            switch(stiffness_type) {
            case StiffnessType::TRIAL:
                D->assemble_trial_geometry();
                break;
            case StiffnessType::CURRENT:
                D->assemble_current_geometry();
                break;
            case StiffnessType::INITIAL:
                // initial geometry mostly like to be empty but may contain something if initially loaded
                D->assemble_initial_geometry();
                break;
            }

            // swap back
            std::swap(t_fox, t_geometry);

            return t_fox;
        });

        // handle stiffness matrix
        // backup whatever
        std::swap(access::rw(current_stiffness), t_stiff);

        // now t_stiffness is empty
        t_stiff = triplet_form<double, uword>(n_block, n_block, num_entry);

        switch(stiffness_type) {
        case StiffnessType::TRIAL:
            D->assemble_trial_stiffness();
            break;
        case StiffnessType::CURRENT:
            D->assemble_current_stiffness();
            break;
        case StiffnessType::INITIAL:
            D->assemble_initial_stiffness();
            break;
        }

        // account for geometry matrix
        t_stiff += f_geometry.get();

        // need to be performed before applying constraints
        f_mass.get();

        // now apply constraints
        if(SUANPAN_SUCCESS != LeeNewmarkBase::process_constraint()) return SUANPAN_FAIL;
        t_stiff.csc_condense();

        // move original stiffness matrix back
        // need to add it to global stiffness matrix later
        std::swap(access::rw(current_stiffness), t_stiff);

        // now current mass and stiffness are formulated
        // assemble unrolled damping matrix and the corresponding damping force
        update_stiffness();
        update_residual();

        if(StiffnessType::TRIAL != stiffness_type) {
            const auto& row = t_triplet.row_mem();
            const auto& col = t_triplet.col_mem();
            const auto& val = t_triplet.val_mem();

            // left top block of unrolled damping matrix may not be zero
            // make a copy
            auto& t_rabbit = access::rw(rabbit);
            t_rabbit = triplet_form<double, uword>(n_block, n_block, t_stiff.n_elem);
            for(uword I = 0; I < t_triplet.n_elem; ++I) {
                // quit if current column is beyond the original size of matrix
                if(col[I] >= n_block) break;
                // check in left top block of unrolled damping matrix to be used in subsequent iterations
                if(row[I] < n_block) t_rabbit.at(row[I], col[I]) = val[I];
            }
        }

        stiffness += t_stiff;

        first_iteration = false;
    }
    else {
        // if not first iteration
        // erase the tangent stiffness entries

        uword *ptr_a, *ptr_b;

        if(t_triplet.is_csc_sorted()) {
            ptr_a = t_triplet.col_mem();
            ptr_b = t_triplet.row_mem();
        }
        else if(t_triplet.is_csr_sorted()) {
            ptr_a = t_triplet.row_mem();
            ptr_b = t_triplet.col_mem();
        }
        else {
            suanpan_error("the system is not sorted while entering iteration, please file a bug report.\n");
            return SUANPAN_FAIL;
        }

        const auto& val = t_triplet.val_mem();

        for(uword I = 0; I < t_triplet.n_elem; ++I) {
            // quit if current column/row is beyond the original size of matrix
            if(ptr_a[I] >= n_block) break;
            // erase existing entries if fall in intact stiffness matrix
            if(ptr_b[I] < n_block) val[I] = 0.;
        }

        // check in original nonzero entries in unrolled damping matrix
        stiffness += rabbit;

        update_residual();

        stiffness += t_stiff;
    }

    return SUANPAN_SUCCESS;
}

void LeeNewmarkFull::print() { suanpan_info("A Newmark solver using Lee's damping model with adjustable bandwidth using %s stiffness. doi: 10.1016/j.compstruc.2020.106423 and 10.1016/j.compstruc.2021.106663\n", stiffness_type == StiffnessType::TRIAL ? "tangent" : stiffness_type == StiffnessType::CURRENT ? "converged" : "initial"); }
