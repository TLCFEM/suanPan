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

#include "ModalRecorder.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <magic_enum/magic_enum.hpp>
#ifdef SUANPAN_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

extern fs::path SUANPAN_OUTPUT;

void ModalRecorder::record_impl(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    frequency.reset();
    participation_factor.reset();
    effective_mass.reset();
    mode_shape.reset();

    node_dof.clear();

    if(W->get_analysis_type() != AnalysisType::EIGEN) return;

    auto& eigenvalue = W->get_eigenvalue();
    auto& eigenvector = W->get_eigenvector();
    if(eigenvalue.is_empty() || eigenvector.is_empty()) return;

    const auto mass = W->get_mass();
    if(nullptr == mass || mass->is_empty()) {
        suanpan_warning("The mass matrix is not available, modal properties cannot be computed by recorder {}.\n", get_tag());
        return;
    }

    const auto n_mode = std::min(eigenvalue.n_elem, eigenvector.n_cols);
    const auto n_dof = eigenvector.n_rows;
    if(0 == n_mode || 0 == n_dof) return;

    const auto components = static_cast<uword>(directions.size());

    // unit influence vectors in all requested directions
    mat influence(n_dof, components, fill::zeros);

    for(auto& I : D->get_node_pool()) {
        if(!I->is_active()) continue;

        auto& dof = I->get_reordered_dof();
        if(dof.is_empty()) continue;

        node_dof.emplace_back(I->get_tag(), dof);

        for(auto J{0u}; J < components; ++J)
            for(const auto K : I->get_dof({directions[J]})) influence(K, J) = 1.;
    }

    // M * phi and M * r, computed once for all modes and directions
    const mat mass_eigenvector = *mass * (mode_shape = eigenvector.cols(0, n_mode - 1));
    const mat mass_influence = *mass * influence;

    frequency.zeros(n_mode);
    participation_factor.zeros(n_mode, components);
    effective_mass.zeros(n_mode, components);

    const auto mass_threshold = 1E-12 * std::max(1., mass->max());

    for(uword I{0}; I < n_mode; ++I) {
        frequency(I) = std::sqrt(std::max(0., eigenvalue(I))) / datum::tau;

        const auto modal_mass = dot(mode_shape.col(I), mass_eigenvector.col(I));

        // massless modes carry no modal information
        if(modal_mass <= mass_threshold) continue;

        // normalize the mode shape with respect to the mass matrix to make the participation factor unique
        mode_shape.col(I) /= std::sqrt(modal_mass);

        for(auto J{0u}; J < components; ++J) {
            const auto gamma = dot(mode_shape.col(I), mass_influence.col(J));
            participation_factor(I, J) = gamma;
            effective_mass(I, J) = gamma * gamma;
        }
    }
}

ModalRecorder::ModalRecorder(const unsigned T, std::vector<Node::DOF>&& D, const bool H)
    : Recorder(T, {}, OutputType::NL, 1, H)
    , directions(std::move(D)) {
    std::ranges::sort(directions);
    const auto [first, last] = std::ranges::unique(directions);
    directions.erase(first, last);
}

void ModalRecorder::save() {
    if(frequency.is_empty()) return;

#ifdef SUANPAN_HDF5
    std::ostringstream file_name;
    if(comm_size > 1) file_name << 'P' << comm_rank << '-';
    file_name << 'R' << get_tag() << "-Modal.h5";

    const auto file_id = H5Fcreate((SUANPAN_OUTPUT / file_name.str()).generic_string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimension[2]{frequency.n_elem, 1};
    H5LTmake_dataset(file_id, "Frequency", 2, dimension, H5T_NATIVE_DOUBLE, frequency.mem);

    for(std::size_t I{0}; I < directions.size(); ++I) {
        const auto group_name = std::string(magic_enum::enum_name(directions[I]));
        const auto group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        dimension[0] = participation_factor.n_rows;
        H5LTmake_dataset(group_id, "ParticipationFactor", 2, dimension, H5T_NATIVE_DOUBLE, participation_factor.colptr(I));
        H5LTmake_dataset(group_id, "EffectiveMass", 2, dimension, H5T_NATIVE_DOUBLE, effective_mass.colptr(I));

        H5Gclose(group_id);
    }

    for(uword I{0}; I < mode_shape.n_cols; ++I) {
        const auto group_id = H5Gcreate(file_id, ("Mode " + std::to_string(I + 1)).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        const vec mode_shape_column = mode_shape.col(I);

        for(const auto& [tag, dof] : node_dof) {
            const vec nodal_shape = mode_shape_column(dof);
            dimension[0] = nodal_shape.n_elem;
            dimension[1] = 1;
            H5LTmake_dataset(group_id, ("N" + std::to_string(tag)).c_str(), 2, dimension, H5T_NATIVE_DOUBLE, nodal_shape.mem);
        }

        H5Gclose(group_id);
    }

    H5Fclose(file_id);
#endif
}

void ModalRecorder::print() { suanpan_info("A recorder to record modal properties, including frequencies, mode shapes, modal participation factors and effective modal masses.\n"); }
