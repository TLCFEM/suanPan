/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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
#include <Step/Frequency.h>

/**
 * \brief method to apply the BC to the system.
 * \param D `Domain`
 * \return 0
 */
int MultiplierBC::process(const shared_ptr<DomainBase>& D) {
	if(auto& t_step = *D->get_current_step(); typeid(Frequency) == typeid(t_step)) return PenaltyBC::process(D);

	auto& t_stiff = D->get_factory()->get_stiffness();

	for(const auto I : get_nodal_active_dof(D)) {
		D->insert_restrained_dof(I);
		t_stiff->unify(I);
	}

	return SUANPAN_SUCCESS;
}
