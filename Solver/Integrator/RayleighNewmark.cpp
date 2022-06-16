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

#include "RayleighNewmark.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Element/Utility/MatrixModifier.hpp>

RayleighNewmark::RayleighNewmark(const unsigned T, const double A, const double B, const double DA, const double DB, const double DC, const double DD)
    : Newmark(T, A, B)
    , damping_alpha(DA)
    , damping_beta(DB)
    , damping_zeta(DC)
    , damping_eta(DD) {}

void RayleighNewmark::assemble_resistance() {
    suanpan::for_all(get_domain().lock()->get_element_pool(), [&](const shared_ptr<Element>& t_element) { suanpan::damping::rayleigh::apply(t_element, damping_alpha, damping_beta, damping_zeta, damping_eta); });

    Newmark::assemble_resistance();
}
