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

#include "EigenRecorder.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Toolbox/HDF5Wrapper.h>

extern fs::path SUANPAN_OUTPUT;

void EigenRecorder::record_impl(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    if(W->get_analysis_type() != AnalysisType::EIGEN) return;

    eigen_value = W->get_eigenvalue();
    const auto& eig_vector = W->get_eigenvector();

    eigen_pool = std::vector<std::map<unsigned, vec>>(eigen_value.n_elem);

    const auto& node_pool = D->get_node_pool();

    for(uword I{0}; I < eigen_value.n_elem; ++I) {
        const vec t_eigen_vector = eig_vector.col(I);
        for(const auto& J : node_pool) eigen_pool[I].emplace(J->get_tag(), t_eigen_vector(J->get_reordered_dof()));
    }
}

EigenRecorder::EigenRecorder(const unsigned T, const bool H)
    : Recorder(T, {}, OutputType::NL, 1, H) {}

void EigenRecorder::save(const shared_ptr<DomainBase>&) {
    if(eigen_value.is_empty()) return;

#ifdef SUANPAN_HDF5
    std::ostringstream file_name;
    // ReSharper disable once CppIfCanBeReplacedByConstexprIf
    // ReSharper disable once CppUnreachableCode
    if(comm_size > 1) file_name << 'P' << comm_rank << '-';
    file_name << "Eigenvalue.h5";

    hdf5::File file(SUANPAN_OUTPUT / file_name.str(), hdf5::File::Mode::Truncate);

    file.create_dataset("Eigenvalue", eigen_value.mem, {static_cast<hsize_t>(eigen_value.n_elem), 1});

    for(uword I{0}, J{1}; I < eigen_value.n_elem; ++I, ++J) {
        auto group = file.create_group("Eigenvalue " + std::to_string(J));

        for(const auto& [e_val, e_vec] : eigen_pool[I]) group.create_dataset('N' + std::to_string(e_val), e_vec.mem, {static_cast<hsize_t>(e_vec.n_cols), static_cast<hsize_t>(e_vec.n_rows)});
    }
#endif
}

void EigenRecorder::print() { suanpan_info("A recorder to record eigenvalues and eigenvectors.\n"); }
