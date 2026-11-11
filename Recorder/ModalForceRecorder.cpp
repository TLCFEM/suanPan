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

    mode_count = 0;
    node_force.clear();
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

    auto mode_shape = eigenvector;
    const mat mass_eigenvector = *mass * mode_shape;
    const auto mass_threshold = 1E-12 * std::max(1., mass->max());

    for(uword I{0}; I < n_mode; ++I)
        if(const auto modal_mass = dot(mode_shape.col(I), mass_eigenvector.col(I)); modal_mass <= mass_threshold) mode_shape.col(I).zeros();
        else mode_shape.col(I) /= std::sqrt(modal_mass);

    mat force_response(n_dof, n_mode, fill::zeros);

    element_force.reserve(D->get_element_pool().size() / comm_size + 1);
    for(const auto& t_ele : D->get_element_pool()) {
        if(!t_ele->is_local) continue;

        auto& stiffness = t_ele->get_trial_stiffness();
        if(stiffness.is_empty()) continue;

        auto& encoding = t_ele->get_dof_encoding();
        mat local_force = stiffness * mode_shape.rows(encoding);
        force_response.rows(encoding) += local_force;
        element_force.emplace_back(t_ele->get_tag(), std::move(local_force));
    }

    std::ranges::sort(element_force, [](const auto& A, const auto& B) { return A.first < B.first; });

    allreduce(force_response);

    // store only the nodal quantities that are actually written to disk
    node_force.reserve(D->get_node_pool().size());
    for(auto& I : D->get_node_pool()) node_force.emplace_back(I->get_tag(), force_response.rows(I->get_reordered_dof()));

    mode_count = n_mode;
}

ModalForceRecorder::ModalForceRecorder(const unsigned T, const bool H)
    : Recorder(T, {}, OutputType::NL, 1, H) {}

void ModalForceRecorder::save() {
    if(0 == mode_count) return;

#ifdef SUANPAN_HDF5
    std::ostringstream file_name;
    // ReSharper disable once CppIfCanBeReplacedByConstexprIf
    // ReSharper disable once CppUnreachableCode
    if(comm_size > 1) file_name << 'P' << comm_rank << '-';
    file_name << 'R' << get_tag() << "-ModalForce.h5";

    const auto file_id = H5Fcreate((SUANPAN_OUTPUT / file_name.str()).generic_string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimension[2]{1, 1};

    for(uword I{0}; I < mode_count; ++I) {
        const auto mode_id = H5Gcreate(file_id, ("Mode " + std::to_string(I + 1)).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        const auto nodes_id = H5Gcreate(mode_id, "Nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        for(const auto& [tag, nodal_force] : node_force) {
            dimension[0] = nodal_force.n_rows;

            H5LTmake_dataset(nodes_id, ("N" + std::to_string(tag)).c_str(), 2, dimension, H5T_NATIVE_DOUBLE, nodal_force.colptr(I));
        }
        H5Gclose(nodes_id);

        const auto elements_id = H5Gcreate(mode_id, "Elements", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        for(const auto& [tag, elemental_force] : element_force) {
            dimension[0] = elemental_force.n_rows;

            H5LTmake_dataset(elements_id, ("E" + std::to_string(tag)).c_str(), 2, dimension, H5T_NATIVE_DOUBLE, elemental_force.colptr(I));
        }
        H5Gclose(elements_id);

        H5Gclose(mode_id);
    }

    H5Fclose(file_id);
#endif
}

void ModalForceRecorder::print() { suanpan_info("A recorder to record modal elastic forces, including the assembled nodal force vector and the element force vectors.\n"); }
