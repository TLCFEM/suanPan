/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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

#include "ShellBase.h"
#include <Toolbox/tensor.h>

const uvec ShellBase::m_dof{0, 1, 5};
const uvec ShellBase::p_dof{2, 3, 4};

vec ShellBase::reshuffle(const vec& membrane_resistance, const vec& plate_resistance) {
    suanpan_assert([&] { if(membrane_resistance.n_elem != plate_resistance.n_elem) throw std::invalid_argument("size conflicts"); });

    const auto t_size = 2 * plate_resistance.n_elem;

    vec total_resistance(t_size);

    for(unsigned I = 0, J = 0; I < t_size; I += 6, J += 3) {
        const span K(J, J + 2llu);
        total_resistance(I + m_dof) = membrane_resistance(K);
        total_resistance(I + p_dof) = plate_resistance(K);
    }

    return total_resistance;
}

mat ShellBase::reshuffle(const mat& membrane_stiffness, const mat& plate_stiffness, const mat& mp_stiffness, const mat& pm_stiffness) {
    const auto t_size = 2 * plate_stiffness.n_cols;

    mat total_stiffness(t_size, t_size, fill::zeros);

    for(auto J = 0llu, L = 0llu; J < t_size; J += 6llu, L += 3llu) {
        const span N(L, L + 2llu);
        for(auto I = 0llu, K = 0llu; I < t_size; I += 6llu, K += 3llu) {
            const span M(K, K + 2llu);
            total_stiffness(I + m_dof, J + m_dof) = membrane_stiffness(M, N);
            total_stiffness(I + p_dof, J + p_dof) = plate_stiffness(M, N);
            total_stiffness(I + m_dof, J + p_dof) = mp_stiffness(M, N);
            total_stiffness(I + p_dof, J + m_dof) = pm_stiffness(M, N);
        }
    }

    return total_stiffness;
}

void ShellBase::direction_cosine() {
    if(!is_nlgeom() && !trans_mat.is_empty()) return;

    const mat coor = get_coordinate(3).t();

    trans_mat.set_size(3, 3);

    vec x1 = coor.col(1) - coor.col(0);
    vec x2 = coor.col(2) - coor.col(0);

    if(is_nlgeom()) {
        const vec disp = get_trial_displacement();
        x1 += disp.rows(6, 8) - disp.rows(0, 2);
        x2 += disp.rows(12, 14) - disp.rows(0, 2);
    }

    trans_mat.col(0) = normalise(x1);
    trans_mat.col(2) = normalise(cross(x1, x2));
    trans_mat.col(1) = cross(trans_mat.col(2), trans_mat.col(0));
}

mat ShellBase::get_local_coordinate() const {
    // mat l_coordinate = get_coordinate(3).t();
    // l_coordinate = trans_mat.t() * (l_coordinate - repmat(l_coordinate.col(0), 1, l_coordinate.n_cols));
    // if(norm(l_coordinate.row(2)) > 1E-6) suanpan_warning("non-planar shell geometry detected.\n");
    // return l_coordinate.rows(0, 1).t();
    return (get_coordinate(3) * trans_mat).eval().cols(0, 1);
}

vec& ShellBase::transform_from_local_to_global(vec& resistance) const {
    for(auto I = 0llu, J = 2llu; I < resistance.n_elem; I += 3llu, J += 3llu) resistance.rows(I, J) = trans_mat * resistance.rows(I, J);

    return resistance;
}

vec& ShellBase::transform_from_global_to_local(vec& displacement) const {
    for(auto I = 0llu, J = 2llu; I < displacement.n_elem; I += 3llu, J += 3llu) displacement.rows(I, J) = trans_mat.t() * displacement.rows(I, J);

    return displacement;
}

mat& ShellBase::transform_from_local_to_global(mat& stiffness) const {
    suanpan_assert([&] { if(stiffness.n_cols != stiffness.n_rows) throw std::invalid_argument("size conflicts"); });

    mat global_trans(size(stiffness), fill::zeros);
    for(auto I = 0llu, J = 2llu; I < stiffness.n_cols; I += 3llu, J += 3llu) {
        const span t_span(I, J);
        global_trans(t_span, t_span) = trans_mat;
    }

    return stiffness = global_trans * stiffness * global_trans.t();
}

vec ShellBase::transform_from_global_to_local(const vec& displacement) const {
    auto transformed = displacement;
    transform_from_global_to_local(transformed);
    return transformed;
}

mat ShellBase::transform_from_local_to_global(mat&& stiffness) const { return std::move(transform_from_local_to_global(stiffness)); }

mat ShellBase::transform_to_global_geometry(const mat& stiffness, const vec& resistance, const vec& displacement) const {
    static const span a(0, 2), b(3, 5), c(6, 8);
    static const uvec d{0, 1, 2, 6, 7, 8, 12, 13, 14};

    const mat coor = get_coordinate(3).t();
    const vec3 x1 = coor.col(0) + displacement(span(0, 2));
    const vec3 x2 = coor.col(1) + displacement(span(6, 8));
    const vec3 x3 = coor.col(2) + displacement(span(12, 14));
    const auto diff_triad = tensor::diff_triad(x1, x2, x3);

    mat left(size(stiffness), fill::zeros), right(size(stiffness), fill::zeros);

    for(auto I = 0llu, J = 1llu, K = 2llu; I < resistance.n_elem; I += 3llu, J += 3llu, K += 3llu) {
        left(uvec{I, J, K}, d) = resistance(I) * diff_triad.rows(a) + resistance(J) * diff_triad.rows(b) + resistance(K) * diff_triad.rows(c);
        const vec nodal_disp = displacement(span(I, K));
        right(d, uvec{I}) = diff_triad.rows(a).t() * nodal_disp;
        right(d, uvec{J}) = diff_triad.rows(b).t() * nodal_disp;
        right(d, uvec{K}) = diff_triad.rows(c).t() * nodal_disp;
    }

    mat global_trans(size(stiffness), fill::zeros);
    for(auto I = 0llu, J = 2llu; I < stiffness.n_cols; I += 3llu, J += 3llu) {
        const span t_span(I, J);
        global_trans(t_span, t_span) = trans_mat;
    }

    return left + global_trans * stiffness * right.t();
}
