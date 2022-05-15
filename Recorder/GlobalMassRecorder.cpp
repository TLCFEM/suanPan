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

#include "GlobalMassRecorder.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Element/Element.h>

void GlobalMassRecorder::assemble_mass(const mat& EM, const uvec& EC, mat& GM) {
    if(EM.is_empty()) return;

    for(unsigned I = 0; I < EC.n_elem; ++I) for(unsigned J = 0; J < EC.n_elem; ++J) GM(EC(J), EC(I)) += EM(J, I);
}

GlobalMassRecorder::GlobalMassRecorder(const unsigned T, const unsigned I, const bool R, const bool H)
    : GlobalRecorder(T, OutputType::M, I, R, H) {}

void GlobalMassRecorder::record(const shared_ptr<DomainBase>& D) {
    if(1 != interval && counter++ != interval) return;

    counter = 1;

    auto& W = D->get_factory();
    auto& C = D->get_color_map();

    const uword S = W->get_size();

    vec mass(S * S, fill::zeros);

    if(mat g_mass(mass.memptr(), S, S, false, true); C.empty()) for(const auto& I : D->get_element_pool()) assemble_mass(I->get_current_mass(), I->get_dof_encoding(), g_mass);
    else
        std::ranges::for_each(C, [&](const vector<unsigned>& color) {
            suanpan_for_each(color.begin(), color.end(), [&](const unsigned tag) {
                const auto& I = D->get<Element>(tag);
                assemble_mass(I->get_current_mass(), I->get_dof_encoding(), g_mass);
            });
        });

    insert({mass}, 0);

    if(if_record_time()) insert(D->get_factory()->get_current_time());
}

void GlobalMassRecorder::print() { suanpan_info("A Global Mass Recorder.\n"); }
