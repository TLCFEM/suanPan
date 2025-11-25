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

#include "GlobalRecorder.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Element/Element.h>

void GlobalRecorder::assemble_matrix(const mat& local, const uvec& encoding, mat& global) {
    if(local.is_empty()) return;

    for(unsigned I = 0; I < encoding.n_elem; ++I)
        for(unsigned J = 0; J < encoding.n_elem; ++J) global(encoding(J), encoding(I)) += local(J, I);
}

GlobalRecorder::GlobalRecorder(const unsigned T, const OutputType L, const unsigned I, const bool H)
    : Recorder(T, {0}, L, I, H) {}

void GlobalRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    if(OutputType::KE == variable_type) {
        auto kinetic_energy = 0.;
        for(auto& I : D->get_pool<Element>()) kinetic_energy += I->get_kinetic_energy();
        insert({{kinetic_energy, D->get_factory()->get_kinetic_energy()}}, 0);
    }
    else if(OutputType::VE == variable_type) {
        auto viscous_energy = 0.;
        for(auto& I : D->get_pool<Element>()) viscous_energy += I->get_viscous_energy();
        insert({{viscous_energy, D->get_factory()->get_viscous_energy()}}, 0);
    }
    else if(OutputType::NVE == variable_type) {
        auto nonviscous_energy = 0.;
        for(auto& I : D->get_pool<Element>()) nonviscous_energy += I->get_nonviscous_energy();
        insert({{nonviscous_energy, D->get_factory()->get_nonviscous_energy()}}, 0);
    }
    else if(OutputType::SE == variable_type) {
        auto strain_energy = 0.;
        for(auto& I : D->get_pool<Element>()) strain_energy += I->get_strain_energy();
        insert({{strain_energy, D->get_factory()->get_strain_energy()}}, 0);
    }
    else insert({{.0, .0}}, 0);

    insert(D->get_factory()->get_current_time());
}

void GlobalRecorder::print() { suanpan_info("A global recorder.\n"); }

GlobalStiffnessRecorder::GlobalStiffnessRecorder(const unsigned T, const unsigned I, const bool H)
    : GlobalRecorder(T, OutputType::K, I, H) {}

void GlobalStiffnessRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    auto& W = D->get_factory();
    auto& C = D->get_color_map();

    const uword S = W->get_size();

    vec stiffness(S * S, fill::zeros);

    if(mat g_stiffness(stiffness.memptr(), S, S, false, true); C.empty())
        for(const auto& I : D->get_element_pool()) assemble_matrix(I->get_current_stiffness(), I->get_dof_encoding(), g_stiffness);
    else
        std::ranges::for_each(C, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = D->get<Element>(tag);
                assemble_matrix(I->get_current_stiffness(), I->get_dof_encoding(), g_stiffness);
            });
        });

    insert({stiffness}, 0);

    insert(D->get_factory()->get_current_time());
}

void GlobalStiffnessRecorder::print() { suanpan_info("A global stiffness recorder.\n"); }

GlobalMassRecorder::GlobalMassRecorder(const unsigned T, const unsigned I, const bool H)
    : GlobalRecorder(T, OutputType::M, I, H) {}

void GlobalMassRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    auto& W = D->get_factory();
    auto& C = D->get_color_map();

    const uword S = W->get_size();

    vec mass(S * S, fill::zeros);

    if(mat g_mass(mass.memptr(), S, S, false, true); C.empty())
        for(const auto& I : D->get_element_pool()) assemble_matrix(I->get_current_mass(), I->get_dof_encoding(), g_mass);
    else
        std::ranges::for_each(C, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = D->get<Element>(tag);
                assemble_matrix(I->get_current_mass(), I->get_dof_encoding(), g_mass);
            });
        });

    insert({mass}, 0);

    insert(D->get_factory()->get_current_time());
}

void GlobalMassRecorder::print() { suanpan_info("A global mass recorder.\n"); }
