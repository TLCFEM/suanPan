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

#include <Domain/DOF.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Element/Element.h>

GlobalRecorder::GlobalRecorder(const unsigned T, const OutputType L, const unsigned I, const bool R, const bool H)
    : Recorder(T, {0}, L, I, R, H) {}

void GlobalRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    auto get_momentum_component = [&](const DOF C) {
        auto momentum = 0.;
        for(auto& I : D->get_pool<Element>()) momentum += I->get_momentum_component(C);
        return momentum;
    };

    if(OutputType::KE == get_variable_type()) {
        auto kinetic_energy = 0.;
        for(auto& I : D->get_pool<Element>()) kinetic_energy += I->get_kinetic_energy();
        insert({{kinetic_energy, D->get_factory()->get_kinetic_energy()}}, 0);
    }
    else if(OutputType::VE == get_variable_type()) {
        auto viscous_energy = 0.;
        for(auto& I : D->get_pool<Element>()) viscous_energy += I->get_viscous_energy();
        insert({{viscous_energy, D->get_factory()->get_viscous_energy()}}, 0);
    }
    else if(OutputType::NVE == get_variable_type()) {
        auto nonviscous_energy = 0.;
        for(auto& I : D->get_pool<Element>()) nonviscous_energy += I->get_nonviscous_energy();
        insert({{nonviscous_energy, D->get_factory()->get_nonviscous_energy()}}, 0);
    }
    else if(OutputType::SE == get_variable_type()) {
        auto strain_energy = 0.;
        for(auto& I : D->get_pool<Element>()) strain_energy += I->get_strain_energy();
        insert({{strain_energy, D->get_factory()->get_strain_energy()}}, 0);
    }
    else if(OutputType::MM == get_variable_type()) {
        auto momentum = 0.;
        for(auto& I : D->get_pool<Element>()) momentum += accu(I->get_momentum());
        insert({{momentum, accu(D->get_factory()->get_momentum())}}, 0);
    }
    else if(OutputType::MM1 == get_variable_type()) insert({{get_momentum_component(DOF::U1)}}, 0);
    else if(OutputType::MM2 == get_variable_type()) insert({{get_momentum_component(DOF::U2)}}, 0);
    else if(OutputType::MM3 == get_variable_type()) insert({{get_momentum_component(DOF::U3)}}, 0);
    else if(OutputType::MM4 == get_variable_type() || OutputType::MMR1 == get_variable_type()) insert({{get_momentum_component(DOF::UR1)}}, 0);
    else if(OutputType::MM5 == get_variable_type() || OutputType::MMR2 == get_variable_type()) insert({{get_momentum_component(DOF::UR2)}}, 0);
    else if(OutputType::MM6 == get_variable_type() || OutputType::MMR3 == get_variable_type()) insert({{get_momentum_component(DOF::UR3)}}, 0);
    else insert({{.0, .0}}, 0);

    if(if_record_time()) insert(D->get_factory()->get_current_time());
}

void GlobalRecorder::print() {
    suanpan_info("A global recorder.\n");
}
