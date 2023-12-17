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

#include "MultiplierBC.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Step/Step.h>
#include <Solver/Integrator/Integrator.h>

/**
 * \brief method to apply the BC to the system.
 * \param D `Domain`
 * \return 0
 */
int MultiplierBC::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    // record restrained DoFs to erase machine error
    // the container used is concurrently safe
    D->insert_restrained_dof(dof_encoding = get_nodal_active_dof(D));

    if(IntegratorType::Explicit == D->get_current_step()->get_integrator()->type()) {
        if(auto& t_mass = W->get_mass(); nullptr != t_mass) {
            std::scoped_lock lock(W->get_mass_mutex());
            for(const auto I : dof_encoding) t_mass->unify(I);
        }
    }
    else {
        if(auto& t_stiff = W->get_stiffness(); nullptr != t_stiff) {
            std::scoped_lock lock(W->get_stiffness_mutex());
            for(const auto I : dof_encoding) t_stiff->unify(I);
        }
        if(auto& t_mass = W->get_mass(); nullptr != t_mass) {
            std::scoped_lock lock(W->get_mass_mutex());
            for(const auto I : dof_encoding) t_mass->nullify(I);
        }
        if(auto& t_damping = W->get_damping(); nullptr != t_damping) {
            std::scoped_lock lock(W->get_damping_mutex());
            for(const auto I : dof_encoding) t_damping->nullify(I);
        }
        if(auto& t_nonviscous = W->get_nonviscous(); nullptr != t_nonviscous) {
            std::scoped_lock lock(W->get_nonviscous_mutex());
            for(const auto I : dof_encoding) t_nonviscous->nullify(I);
        }
        if(auto& t_geometry = W->get_geometry(); nullptr != t_geometry) {
            std::scoped_lock lock(W->get_geometry_mutex());
            for(const auto I : dof_encoding) t_geometry->nullify(I);
        }
    }

    return SUANPAN_SUCCESS;
}

int MultiplierBC::process_resistance(const shared_ptr<DomainBase>& D) {
    D->insert_restrained_dof(dof_encoding);

    return SUANPAN_SUCCESS;
}
