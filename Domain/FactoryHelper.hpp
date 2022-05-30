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

#ifndef FACTORY_HELPER_HPP
#define FACTORY_HELPER_HPP

#include <suanPan.h>

#include "Factory.hpp"
#include "MetaMat/MetaMat"

template<sp_d T> unique_ptr<MetaMat<T>> get_basic_container(const Factory<T>* const W) {
    switch(W->storage_type) {
    case StorageScheme::FULL:
#ifdef SUANPAN_CUDA
        if(SolverType::CUDA == W->solver) return std::make_unique<FullMatCUDA<T>>(W->n_size, W->n_size);
#endif
        return std::make_unique<FullMat<T>>(W->n_size, W->n_size);
    case StorageScheme::BAND:
        if(SolverType::SPIKE == W->solver) return std::make_unique<BandMatSpike<T>>(W->n_size, W->n_lobw, W->n_upbw);
        return std::make_unique<BandMat<T>>(W->n_size, W->n_lobw, W->n_upbw);
    case StorageScheme::BANDSYMM:
        return std::make_unique<BandSymmMat<T>>(W->n_size, W->n_lobw);
    case StorageScheme::SYMMPACK:
        return std::make_unique<SymmPackMat<T>>(W->n_size);
    case StorageScheme::SPARSE:
        if(SolverType::MUMPS == W->solver) return std::make_unique<SparseMatMUMPS<T>>(W->n_size, W->n_size, W->n_elem);
        if(SolverType::SUPERLU == W->solver) return std::make_unique<SparseMatSuperLU<T>>(W->n_size, W->n_size, W->n_elem);
#ifdef SUANPAN_MKL
        if(SolverType::PARDISO == W->solver) return std::make_unique<SparseMatPARDISO<T>>(W->n_size, W->n_size, W->n_elem);
        if(SolverType::FGMRES == W->solver) return std::make_unique<SparseMatFGMRES<T>>(W->n_size, W->n_size, W->n_elem);
#endif
#ifdef SUANPAN_CUDA
        if(SolverType::CUDA == W->solver) return std::make_unique<SparseMatCUDA<T>>(W->n_size, W->n_size, W->n_elem);
#endif
        return std::make_unique<SparseMatSuperLU<T>>(W->n_size, W->n_size, W->n_elem);
    case StorageScheme::SPARSESYMM:
#ifdef SUANPAN_MKL
        if(SolverType::FGMRES == W->solver) return std::make_unique<SparseSymmMatFGMRES<T>>(W->n_size, W->n_size, W->n_elem);
#endif
        return std::make_unique<SparseSymmMatMUMPS<T>>(W->n_size, W->n_size, W->n_elem);
    default:
        throw invalid_argument("need a proper storage scheme");
    }
}

template<sp_d T> unique_ptr<MetaMat<T>> get_matrix_container(const Factory<T>* const W) {
    auto global_mat = get_basic_container(W);

    global_mat->set_precision(W->precision);
    global_mat->set_tolerance(W->tolerance);
    global_mat->set_refinement(W->refinement);

    return global_mat;
}

template<sp_d T1> shared_ptr<MetaMat<T1>>& get_mass(const shared_ptr<Factory<T1>>& W) { return W->global_mass; }

template<sp_d T1> shared_ptr<MetaMat<T1>>& get_damping(const shared_ptr<Factory<T1>>& W) { return W->global_damping; }

template<sp_d T1> shared_ptr<MetaMat<T1>>& get_stiffness(const shared_ptr<Factory<T1>>& W) { return W->global_stiffness; }

template<sp_d T1> shared_ptr<MetaMat<T1>>& get_geometry(const shared_ptr<Factory<T1>>& W) { return W->global_geometry; }

template<sp_d T> Col<T>& get_ninja(const shared_ptr<Factory<T>>& W) { return W->ninja; }

template<sp_d T> Col<T>& get_sushi(const shared_ptr<Factory<T>>& W) { return W->sushi; }

template<sp_d T> uvec& get_reference_dof(const shared_ptr<Factory<T>>& W) { return W->reference_dof; }

template<sp_d T> SpMat<T>& get_reference_load(const shared_ptr<Factory<T>>& W) { return W->reference_load; }

template<sp_d T1> uvec& get_auxiliary_encoding(const shared_ptr<Factory<T1>>& W) { return W->auxiliary_encoding; }

template<sp_d T1> Col<T1>& get_auxiliary_lambda(const shared_ptr<Factory<T1>>& W) { return W->auxiliary_lambda; }

template<sp_d T> Col<T>& get_auxiliary_resistance(const shared_ptr<Factory<T>>& W) { return W->auxiliary_resistance; }

template<sp_d T> Col<T>& get_auxiliary_load(const shared_ptr<Factory<T>>& W) { return W->auxiliary_load; }

template<sp_d T> SpMat<T>& get_auxiliary_stiffness(const shared_ptr<Factory<T>>& W) { return W->auxiliary_stiffness; }

template<sp_d T1> SpCol<T1>& get_trial_constraint_resistance(const shared_ptr<Factory<T1>>& W) { return W->trial_constraint_resistance; }

template<sp_d T1> SpCol<T1>& get_current_constraint_resistance(const shared_ptr<Factory<T1>>& W) { return W->current_constraint_resistance; }

template<sp_d T> T& get_trial_time(const shared_ptr<Factory<T>>& W) { return W->trial_time; }

template<sp_d T1> Col<T1>& get_trial_load_factor(const shared_ptr<Factory<T1>>& W) { return W->trial_load_factor; }

template<sp_d T> Col<T>& get_trial_load(const shared_ptr<Factory<T>>& W) { return W->trial_load; }

template<sp_d T> Col<T>& get_trial_settlement(const shared_ptr<Factory<T>>& W) { return W->trial_settlement; }

template<sp_d T> Col<T>& get_trial_resistance(const shared_ptr<Factory<T>>& W) { return W->trial_resistance; }

template<sp_d T> Col<T>& get_trial_damping_force(const shared_ptr<Factory<T>>& W) { return W->trial_damping_force; }

template<sp_d T> Col<T>& get_trial_inertial_force(const shared_ptr<Factory<T>>& W) { return W->trial_inertial_force; }

template<sp_d T> Col<T>& get_trial_displacement(const shared_ptr<Factory<T>>& W) { return W->trial_displacement; }

template<sp_d T> Col<T>& get_trial_velocity(const shared_ptr<Factory<T>>& W) { return W->trial_velocity; }

template<sp_d T> Col<T>& get_trial_acceleration(const shared_ptr<Factory<T>>& W) { return W->trial_acceleration; }

template<sp_d T> Col<T>& get_trial_temperature(const shared_ptr<Factory<T>>& W) { return W->trial_temperature; }

template<sp_d T> T& get_incre_time(const shared_ptr<Factory<T>>& W) { return W->incre_time; }

template<sp_d T1> Col<T1>& get_incre_load_factor(const shared_ptr<Factory<T1>>& W) { return W->incre_load_factor; }

template<sp_d T> Col<T>& get_incre_load(const shared_ptr<Factory<T>>& W) { return W->incre_load; }

template<sp_d T> Col<T>& get_incre_settlement(const shared_ptr<Factory<T>>& W) { return W->incre_settlement; }

template<sp_d T> Col<T>& get_incre_resistance(const shared_ptr<Factory<T>>& W) { return W->incre_resistance; }

template<sp_d T> Col<T>& get_incre_damping_force(const shared_ptr<Factory<T>>& W) { return W->incre_damping_force; }

template<sp_d T> Col<T>& get_incre_inertial_force(const shared_ptr<Factory<T>>& W) { return W->incre_inertial_force; }

template<sp_d T> Col<T>& get_incre_displacement(const shared_ptr<Factory<T>>& W) { return W->incre_displacement; }

template<sp_d T> Col<T>& get_incre_velocity(const shared_ptr<Factory<T>>& W) { return W->incre_velocity; }

template<sp_d T> Col<T>& get_incre_acceleration(const shared_ptr<Factory<T>>& W) { return W->incre_acceleration; }

template<sp_d T> Col<T>& get_incre_temperature(const shared_ptr<Factory<T>>& W) { return W->incre_temperature; }

template<sp_d T> T& get_current_time(const shared_ptr<Factory<T>>& W) { return W->current_time; }

template<sp_d T1> Col<T1>& get_current_load_factor(const shared_ptr<Factory<T1>>& W) { return W->current_load_factor; }

template<sp_d T> Col<T>& get_current_load(const shared_ptr<Factory<T>>& W) { return W->current_load; }

template<sp_d T> Col<T>& get_current_settlement(const shared_ptr<Factory<T>>& W) { return W->current_settlement; }

template<sp_d T> Col<T>& get_current_resistance(const shared_ptr<Factory<T>>& W) { return W->current_resistance; }

template<sp_d T> Col<T>& get_current_damping_force(const shared_ptr<Factory<T>>& W) { return W->current_damping_force; }

template<sp_d T> Col<T>& get_current_inertial_force(const shared_ptr<Factory<T>>& W) { return W->current_inertial_force; }

template<sp_d T> Col<T>& get_current_displacement(const shared_ptr<Factory<T>>& W) { return W->current_displacement; }

template<sp_d T> Col<T>& get_current_velocity(const shared_ptr<Factory<T>>& W) { return W->current_velocity; }

template<sp_d T> Col<T>& get_current_acceleration(const shared_ptr<Factory<T>>& W) { return W->current_acceleration; }

template<sp_d T> Col<T>& get_current_temperature(const shared_ptr<Factory<T>>& W) { return W->current_temperature; }

template<sp_d T> T& get_pre_time(const shared_ptr<Factory<T>>& W) { return W->pre_time; }

template<sp_d T1> Col<T1>& get_pre_load_factor(const shared_ptr<Factory<T1>>& W) { return W->pre_load_factor; }

template<sp_d T> Col<T>& get_pre_load(const shared_ptr<Factory<T>>& W) { return W->pre_load; }

template<sp_d T> Col<T>& get_pre_settlement(const shared_ptr<Factory<T>>& W) { return W->pre_settlement; }

template<sp_d T> Col<T>& get_pre_resistance(const shared_ptr<Factory<T>>& W) { return W->pre_resistance; }

template<sp_d T> Col<T>& get_pre_damping_force(const shared_ptr<Factory<T>>& W) { return W->pre_damping_force; }

template<sp_d T> Col<T>& get_pre_inertial_force(const shared_ptr<Factory<T>>& W) { return W->pre_inertial_force; }

template<sp_d T> Col<T>& get_pre_displacement(const shared_ptr<Factory<T>>& W) { return W->pre_displacement; }

template<sp_d T> Col<T>& get_pre_velocity(const shared_ptr<Factory<T>>& W) { return W->pre_velocity; }

template<sp_d T> Col<T>& get_pre_acceleration(const shared_ptr<Factory<T>>& W) { return W->pre_acceleration; }

template<sp_d T> Col<T>& get_pre_temperature(const shared_ptr<Factory<T>>& W) { return W->pre_temperature; }

template<sp_d T> Col<T>& get_eigenvalue(const shared_ptr<Factory<T>>& W) { return W->eigenvalue; }

template<sp_d T> Mat<T>& get_eigenvector(const shared_ptr<Factory<T>>& W) { return W->eigenvector; }

inline uvec to_uvec(const suanpan_set& in) {
    uvec out(in.size(), fill::none);
    auto I = 0llu;
    for(const auto J : in) out(I++) = J;
    return out;
}

#endif // FACTORY_HELPER_HPP

//! @}
