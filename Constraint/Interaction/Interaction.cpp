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

#include "Interaction.h"

#include <Domain/DomainBase.h>
#include <Element/Element.h>

InteractionPair::InteractionPair(const shared_ptr<Element>& obj_i, const shared_ptr<Element>& obj_j)
    : object_i(obj_i)
    , object_j(obj_j) {}

double InteractionPair::compression() const { return object_j->get(Element::Parameter::RADIUS) + object_i->get(Element::Parameter::RADIUS) - norm(object_j->get_coordinate().t() - object_i->get_coordinate().t() + object_j->get_trial_displacement() - object_i->get_trial_displacement()); }

void Interaction::initialize(const shared_ptr<DomainBase>& D) { domain = D; }
