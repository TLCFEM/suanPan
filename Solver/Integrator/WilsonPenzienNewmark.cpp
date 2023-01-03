/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "WilsonPenzienNewmark.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Toolbox/arpack_wrapper.h>

WilsonPenzienNewmark::WilsonPenzienNewmark(const unsigned T, vec&& DR, const double A, const double B)
    : Newmark(T, A, B)
    , damping_ratio(std::forward<vec>(DR)) {}

int WilsonPenzienNewmark::initialize() {
    if(SUANPAN_SUCCESS != Newmark::initialize()) return SUANPAN_FAIL;

    const auto& D = get_domain();
    auto& W = D->get_factory();

    theta.zeros(W->get_size(), damping_ratio.n_elem);
    beta.zeros(damping_ratio.n_elem);

    return SUANPAN_SUCCESS;
}

int WilsonPenzienNewmark::process_constraint() {
    // process constraint for the first time to obtain proper stiffness
    if(SUANPAN_SUCCESS != Integrator::process_constraint()) return SUANPAN_FAIL;

    auto& W = get_domain()->get_factory();

    auto& t_stiff = W->get_stiffness();
    auto& t_mass = W->get_mass();

    t_stiff->csc_condense();
    t_mass->csc_condense();

    if(first_iteration) {
        cx_vec eig_val;
        cx_mat eig_vec;
        if(SUANPAN_SUCCESS != eig_solve(eig_val, eig_vec, t_stiff->make_copy(), t_mass, static_cast<unsigned>(damping_ratio.n_elem))) {
            if(!eig_pair(eig_val, eig_vec, to_mat(t_stiff), to_mat(t_mass))) {
                suanpan_error("Fail to perform eigen analysis.\n");
                return SUANPAN_FAIL;
            }
            eig_val = eig_val.head(damping_ratio.n_elem);
            eig_vec = eig_vec.head_cols(damping_ratio.n_elem);
        }

        access::rw(theta) = t_mass * mat(abs(eig_vec));
        access::rw(beta) = 2. * damping_ratio % sqrt(abs(eig_val));

        access::rw(first_iteration) = false;
    }

    // the damping matrix is addressed by using the Woodbury formula
    t_stiff += C0 * t_mass + C1 * W->get_damping();

    return SUANPAN_SUCCESS;
}

int WilsonPenzienNewmark::solve(mat& X, const mat& B) {
    mat left, right;

    if(SUANPAN_SUCCESS != Newmark::solve(X, B)) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != Newmark::solve(left, theta)) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != Newmark::solve(right, theta * arma::solve(theta.t() * left + diagmat(1. / C1 / beta), theta.t() * X))) return SUANPAN_FAIL;

    X -= right;

    return SUANPAN_SUCCESS;
}

int WilsonPenzienNewmark::solve(mat& X, const sp_mat& B) {
    mat left, right;

    if(SUANPAN_SUCCESS != Newmark::solve(X, B)) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != Newmark::solve(left, theta)) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != Newmark::solve(right, theta * arma::solve(theta.t() * left + diagmat(1. / C1 / beta), theta.t() * X))) return SUANPAN_FAIL;

    X -= right;

    return SUANPAN_SUCCESS;
}

int WilsonPenzienNewmark::solve(mat& X, mat&& B) {
    mat left, right;

    if(SUANPAN_SUCCESS != Newmark::solve(X, std::forward<mat>(B))) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != Newmark::solve(left, theta)) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != Newmark::solve(right, theta * arma::solve(theta.t() * left + diagmat(1. / C1 / beta), theta.t() * X))) return SUANPAN_FAIL;

    X -= right;

    return SUANPAN_SUCCESS;
}

int WilsonPenzienNewmark::solve(mat& X, sp_mat&& B) {
    mat left, right;

    if(SUANPAN_SUCCESS != Newmark::solve(X, std::forward<sp_mat>(B))) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != Newmark::solve(left, theta)) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != Newmark::solve(right, theta * arma::solve(theta.t() * left + diagmat(1. / C1 / beta), theta.t() * X))) return SUANPAN_FAIL;

    X -= right;

    return SUANPAN_SUCCESS;
}

void WilsonPenzienNewmark::commit_status() {
    first_iteration = true;

    Newmark::commit_status();
}

void WilsonPenzienNewmark::clear_status() {
    first_iteration = true;

    Newmark::clear_status();
}

void WilsonPenzienNewmark::reset_status() {
    first_iteration = true;

    Newmark::reset_status();
}

void WilsonPenzienNewmark::assemble_resistance() {
    const auto& D = get_domain();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_resistance(); });
    auto fb = std::async([&] { D->assemble_damping_force(); }); // consider independent viscous device
    auto fc = std::async([&] { D->assemble_inertial_force(); });

    fa.get();
    fb.get();
    fc.get();

    W->update_trial_damping_force_by(theta * (beta % (theta.t() * W->get_trial_velocity())));

    W->set_sushi(W->get_trial_resistance() + W->get_trial_damping_force() + W->get_trial_inertial_force());
}

void WilsonPenzienNewmark::assemble_matrix() {
    const auto& D = get_domain();
    auto& W = D->get_factory();

    auto fa = std::async([&] { D->assemble_trial_stiffness(); });
    auto fb = std::async([&] { D->assemble_trial_geometry(); });
    auto fc = std::async([&] { D->assemble_trial_damping(); }); // the model may have viscous device
    auto fd = std::async([&] { D->assemble_trial_mass(); });    // need a constant mass matrix

    fa.get();
    fb.get();
    fc.get();
    fd.get();

    auto& t_stiff = W->get_stiffness();

    t_stiff += W->get_geometry();
}

void WilsonPenzienNewmark::print() {
    suanpan_info("A Newmark solver with Wilson-Penzien damping model.\n");
}
