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

#include "MultiplierBC.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

/**
 * \brief method to apply the BC to the system.
 * \param D `Domain`
 * \return 0
 */
int MultiplierBC::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    dof_encoding = get_nodal_active_dof(D);

    for(const auto I : dof_encoding) D->insert_restrained_dof(I);

    if(auto& t_stiff = W->get_stiffness(); nullptr != t_stiff) for(const auto I : dof_encoding) t_stiff->unify(I);
    if(auto& t_mass = W->get_mass(); nullptr != t_mass) for(const auto I : dof_encoding) t_mass->nullify(I);
    if(auto& t_damping = W->get_damping(); nullptr != t_damping) for(const auto I : dof_encoding) t_damping->nullify(I);
    if(auto& t_geometry = W->get_geometry(); nullptr != t_geometry) for(const auto I : dof_encoding) t_geometry->nullify(I);

    return SUANPAN_SUCCESS;
}

int MultiplierBC::process_resistance(const shared_ptr<DomainBase>& D) {
    for(const auto I : get_nodal_active_dof(D)) D->insert_restrained_dof(I);

    return SUANPAN_SUCCESS;
}
