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

#include "Domain.h"

shared_ptr<Amplitude>& get_amplitude(const shared_ptr<Domain>& D, const unsigned T) { return D->amplitude_pond[T]; }

shared_ptr<Expression>& get_expression(const shared_ptr<Domain>& D, const unsigned T) { return D->expression_pond[T]; }

shared_ptr<Constraint>& get_constraint(const shared_ptr<Domain>& D, const unsigned T) { return D->constraint_pond[T]; }

shared_ptr<Converger>& get_converger(const shared_ptr<Domain>& D, const unsigned T) { return D->converger_pond[T]; }

shared_ptr<Criterion>& get_criterion(const shared_ptr<Domain>& D, const unsigned T) { return D->criterion_pond[T]; }

shared_ptr<Database>& get_database(const shared_ptr<Domain>& D, const unsigned T) { return D->database_pond[T]; }

shared_ptr<Element>& get_element(const shared_ptr<Domain>& D, const unsigned T) { return D->element_pond[T]; }

shared_ptr<Group>& get_group(const shared_ptr<Domain>& D, const unsigned T) { return D->group_pond[T]; }

shared_ptr<Integrator>& get_integrator(const shared_ptr<Domain>& D, const unsigned T) { return D->integrator_pond[T]; }

shared_ptr<Load>& get_load(const shared_ptr<Domain>& D, const unsigned T) { return D->load_pond[T]; }

shared_ptr<Material>& get_material(const shared_ptr<Domain>& D, const unsigned T) { return D->material_pond[T]; }

shared_ptr<Modifier>& get_modifier(const shared_ptr<Domain>& D, const unsigned T) { return D->modifier_pond[T]; }

shared_ptr<Node>& get_node(const shared_ptr<Domain>& D, const unsigned T) { return D->node_pond[T]; }

shared_ptr<Orientation>& get_orientation(const shared_ptr<Domain>& D, const unsigned T) { return D->orientation_pond[T]; }

shared_ptr<Recorder>& get_recorder(const shared_ptr<Domain>& D, const unsigned T) { return D->recorder_pond[T]; }

shared_ptr<Section>& get_section(const shared_ptr<Domain>& D, const unsigned T) { return D->section_pond[T]; }

shared_ptr<Solver>& get_solver(const shared_ptr<Domain>& D, const unsigned T) { return D->solver_pond[T]; }

shared_ptr<Step>& get_step(const shared_ptr<Domain>& D, const unsigned T) { return D->step_pond[T]; }

shared_ptr<Amplitude>& get_amplitude(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->amplitude_pond[T]; }

shared_ptr<Expression>& get_expression(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->expression_pond[T]; }

shared_ptr<Constraint>& get_constraint(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->constraint_pond[T]; }

shared_ptr<Converger>& get_converger(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->converger_pond[T]; }

shared_ptr<Criterion>& get_criterion(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->criterion_pond[T]; }

shared_ptr<Database>& get_database(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->database_pond[T]; }

shared_ptr<Element>& get_element(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->element_pond[T]; }

shared_ptr<Group>& get_group(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->group_pond[T]; }

shared_ptr<Integrator>& get_integrator(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->integrator_pond[T]; }

shared_ptr<Load>& get_load(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->load_pond[T]; }

shared_ptr<Material>& get_material(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->material_pond[T]; }

shared_ptr<Modifier>& get_modifier(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->modifier_pond[T]; }

shared_ptr<Node>& get_node(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->node_pond[T]; }

shared_ptr<Orientation>& get_orientation(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->orientation_pond[T]; }

shared_ptr<Recorder>& get_recorder(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->recorder_pond[T]; }

shared_ptr<Section>& get_section(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->section_pond[T]; }

shared_ptr<Solver>& get_solver(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->solver_pond[T]; }

shared_ptr<Step>& get_step(const shared_ptr<DomainBase>& D, const unsigned T) { return std::dynamic_pointer_cast<Domain>(D)->step_pond[T]; }
