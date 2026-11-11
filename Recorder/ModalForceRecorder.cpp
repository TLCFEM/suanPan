/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include "ModalForceRecorder.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Element/Element.h>
#ifdef SUANPAN_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

extern fs::path SUANPAN_OUTPUT;

void ModalForceRecorder::record_impl(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    frequency.reset();
    mode_shape.reset();
    force_response.reset();
    node_dof.clear();
    element_force.clear();

    if(W->get_analysis_type() != AnalysisType::EIGEN) return;

    auto& eigenvalue = W->get_eigenvalue();
    auto& eigenvector = W->get_eigenvector();
    if(eigenvalue.is_empty() || eigenvector.is_empty()) return;

    const auto mass = W->get_mass();
    if(nullptr == mass || mass->is_empty()) {
        suanpan_warning("The mass matrix is not available, modal forces cannot be computed by recorder {}.\n", get_tag());
        return;
    }

    const auto n_mode = std::min(eigenvalue.n_elem, eigenvector.n_cols);
    const auto n_dof = eigenvector.n_rows;
    if(0 == n_mode || 0 == n_dof) return;

    for(auto& I : D->get_node_pool()) {
        if(!I->is_active()) continue;

        auto& dof = I->get_reordered_dof();
        if(dof.is_empty()) continue;

        node_dof.emplace_back(I->get_tag(), dof);
    }

    mode_shape = eigenvector.cols(0, n_mode - 1);

    const mat mass_eigenvector = *mass * mode_shape;
    const auto mass_threshold = 1E-12 * std::max(1., mass->max());

    frequency.zeros(n_mode);

    for(uword I{0}; I < n_mode; ++I) {
        frequency(I) = std::sqrt(std::max(0., eigenvalue(I))) / datum::tau;

        const auto modal_mass = dot(mode_shape.col(I), mass_eigenvector.col(I));

        if(modal_mass <= mass_threshold) {
            mode_shape.col(I).zeros();
            continue;
        }

        mode_shape.col(I) /= std::sqrt(modal_mass);
    }

    force_response.zeros(n_dof, n_mode);

    for(const auto& t_ele : D->get_element_pool()) {
        if(!t_ele->is_local) continue;

        auto& stiffness = t_ele->get_trial_stiffness();
        auto& encoding = t_ele->get_dof_encoding();
        if(stiffness.is_empty() || encoding.n_elem != stiffness.n_rows) continue;

        mat local_force = stiffness * mode_shape.rows(encoding);
        force_response.rows(encoding) += local_force;
        element_force.emplace_back(t_ele->get_tag(), std::move(local_force));
    }

    std::ranges::sort(element_force, [](const auto& A, const auto& B) { return A.first < B.first; });

    allreduce(force_response);
}

ModalForceRecorder::ModalForceRecorder(const unsigned T, const bool H)
    : Recorder(T, {}, OutputType::NL, 1, H) {}

void ModalForceRecorder::save() {
    if(force_response.is_empty()) return;

#ifdef SUANPAN_HDF5
    std::ostringstream file_name;
    // ReSharper disable once CppIfCanBeReplacedByConstexprIf
    // ReSharper disable once CppUnreachableCode
    if(comm_size > 1) file_name << 'P' << comm_rank << '-';
    file_name << 'R' << get_tag() << "-ModalForce.h5";

    const auto file_id = H5Fcreate((SUANPAN_OUTPUT / file_name.str()).generic_string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimension[2]{frequency.n_elem, 1};
    H5LTmake_dataset(file_id, "Frequency", 2, dimension, H5T_NATIVE_DOUBLE, frequency.mem);

    for(uword I{0}; I < mode_shape.n_cols; ++I) {
        const auto group_id = H5Gcreate(file_id, ("Mode " + std::to_string(I + 1)).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for(const auto& [tag, dof] : node_dof) {
            dimension[0] = dof.n_elem;
            dimension[1] = 1;

            const vec nodal_force = force_response(dof, uvec{I});
            H5LTmake_dataset(group_id, ("F" + std::to_string(tag)).c_str(), 2, dimension, H5T_NATIVE_DOUBLE, nodal_force.memptr());
        }

        for(const auto& [tag, e_force] : element_force) {
            dimension[0] = e_force.n_rows;
            dimension[1] = 1;

            const vec elemental_force = e_force.col(I);
            H5LTmake_dataset(group_id, ("EF" + std::to_string(tag)).c_str(), 2, dimension, H5T_NATIVE_DOUBLE, elemental_force.memptr());
        }

        H5Gclose(group_id);
    }

    H5Fclose(file_id);
#endif
}

void ModalForceRecorder::print() { suanpan_info("A recorder to record modal elastic forces, including the assembled nodal force vector and the element force vectors.\n"); }
