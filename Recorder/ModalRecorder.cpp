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
#include <Toolbox/HDF5Wrapper.h>
#endif

extern fs::path SUANPAN_OUTPUT;

void ModalRecorder::record_impl(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    frequency.reset();
    participation_factor.reset();
    effective_mass.reset();
    node_shape.clear();

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

    mat influence(n_dof, components, fill::zeros);
    for(auto& I : D->get_node_pool()) {
        if(!I->is_active()) continue;

        for(auto J{0u}; J < components; ++J)
            for(const auto K : I->get_dof({directions[J]})) influence(K, J) = 1.;
    }

    auto mode_shape = eigenvector;
    const mat mass_eigenvector = *mass * mode_shape;
    const mat mass_influence = *mass * influence;

    frequency.zeros(n_mode);
    participation_factor.zeros(n_mode, components);
    effective_mass.zeros(n_mode, components);

    const auto mass_threshold = 1E-12 * std::max(1., mass->max());
    for(uword I{0}; I < n_mode; ++I) {
        frequency(I) = std::sqrt(std::max(0., eigenvalue(I))) / datum::tau;

        const auto modal_mass = dot(mode_shape.col(I), mass_eigenvector.col(I));
        if(modal_mass <= mass_threshold) {
            mode_shape.col(I).zeros();
            continue;
        }
        mode_shape.col(I) /= std::sqrt(modal_mass);

        for(auto J{0u}; J < components; ++J) {
            const auto gamma = dot(mode_shape.col(I), mass_influence.col(J));
            participation_factor(I, J) = gamma;
            effective_mass(I, J) = gamma * gamma;
        }
    }

    // store only the nodal quantities that are actually written to disk
    node_shape.reserve(D->get_node_pool().size());
    for(auto& I : D->get_node_pool()) node_shape.emplace_back(I->get_tag(), mode_shape.rows(I->get_reordered_dof()));
}

ModalRecorder::ModalRecorder(const unsigned T, std::vector<Node::DOF>&& D, const bool H)
    : Recorder(T, {}, OutputType::NL, 1, H)
    , directions(std::move(D)) {
    std::ranges::sort(directions);
    const auto [first, last] = std::ranges::unique(directions);
    directions.erase(first, last);
}

void ModalRecorder::save(const shared_ptr<DomainBase>&) {
    if(frequency.is_empty()) return;

#ifdef SUANPAN_HDF5
    std::ostringstream file_name;
    // ReSharper disable once CppIfCanBeReplacedByConstexprIf
    // ReSharper disable once CppUnreachableCode
    if(comm_size > 1) file_name << 'P' << comm_rank << '-';
    file_name << 'R' << get_tag() << "-Modal.h5";

    hdf5::File file(SUANPAN_OUTPUT / file_name.str(), hdf5::File::Mode::Truncate);

    file.create_dataset("Frequency", frequency.mem, {static_cast<hsize_t>(frequency.n_elem), 1});

    for(uword I{0}; I < static_cast<uword>(directions.size()); ++I) {
        auto group = file.create_group(magic_enum::enum_name(directions[I]).data());

        const auto rows = static_cast<hsize_t>(participation_factor.n_rows);
        group.create_dataset("ParticipationFactor", participation_factor.colptr(I), {rows, 1});
        group.create_dataset("EffectiveMass", effective_mass.colptr(I), {rows, 1});
    }

    for(uword I{0}; I < frequency.n_elem; ++I) {
        auto group = file.create_group("Mode " + std::to_string(I + 1));

        for(const auto& [tag, shape] : node_shape)
            group.create_dataset("N" + std::to_string(tag), shape.colptr(I), {static_cast<hsize_t>(shape.n_rows), 1});
    }
#endif
}

void ModalRecorder::print() { suanpan_info("A recorder to record modal properties, including frequencies, mode shapes, modal participation factors and effective modal masses.\n"); }
