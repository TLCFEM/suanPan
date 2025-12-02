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
#include <Domain/Factory.hpp>
#include <Element/Element.h>

InteractionPair::InteractionPair(const shared_ptr<Element>& obj_i, const shared_ptr<Element>& obj_j)
    : object_i(obj_i)
    , object_j(obj_j) {
    effective_mass = 1. / (1. / obj_i->get(Element::Parameter::MASS) + 1. / obj_j->get(Element::Parameter::MASS));
    effective_radius = 1. / (1. / obj_i->get(Element::Parameter::RADIUS) + 1. / obj_j->get(Element::Parameter::RADIUS));
    effective_modulus = 1. / ((1. - std::pow(obj_i->get(Element::Parameter::POISSON), 2.)) / obj_i->get(Element::Parameter::ELASTIC) + (1. - std::pow(obj_j->get(Element::Parameter::POISSON), 2.)) / obj_j->get(Element::Parameter::ELASTIC));
    effective_damping = .5 * (obj_i->get(Element::Parameter::DAMPING) + obj_j->get(Element::Parameter::DAMPING));
}

vec InteractionPair::position_i() const { return object_i->get_coordinate(dimension).t() + object_i->get_trial_displacement().head(dimension); }

vec InteractionPair::position_j() const { return object_j->get_coordinate(dimension).t() + object_j->get_trial_displacement().head(dimension); }

vec InteractionPair::velocity_i() const { return object_i->get_trial_velocity().head(dimension); }

vec InteractionPair::velocity_j() const { return object_j->get_trial_velocity().head(dimension); }

const uvec& InteractionPair::dof_i() const { return object_i->get_dof_encoding(); }

const uvec& InteractionPair::dof_j() const { return object_j->get_dof_encoding(); }

double InteractionPair::initial_gap() const { return object_i->get(Element::Parameter::RADIUS) + object_j->get(Element::Parameter::RADIUS); }

void Interaction::initialize(const shared_ptr<DomainBase>& D) { factory = D->get_factory(); }

void Hertzian::apply(const shared_ptr<InteractionPair>& pair) const {
    const vec position_i = pair->position_i(), position_j = pair->position_j();
    const vec chord = position_j - position_i;
    const auto chord_length = norm(chord);
    const auto compression = pair->initial_gap() - chord_length;

    if(compression <= 0.) return;

    const auto normal_factor = 2. * std::sqrt(pair->effective_radius) * pair->effective_modulus * std::pow(compression, .5);
    const auto normal_force_over_length = compression * normal_factor * two_third / chord_length;

    const uvec &dof_i = pair->dof_i(), &dof_j = pair->dof_j();

    {
        const vec repulsive = normal_force_over_length * chord;
        auto& t_resistance = factory->modify_trial_constraint_resistance();
        std::scoped_lock resistance_lock(factory->get_trial_constraint_resistance_mutex());
        for(auto I = 0llu; I < chord.n_elem; ++I) {
            t_resistance(dof_i(I)) += repulsive(I);
            t_resistance(dof_j(I)) -= repulsive(I);
        }
    }

    {
        mat der_repulsive = (normal_factor + normal_force_over_length) / chord_length / chord_length * chord * chord.t();
        der_repulsive.diag() -= normal_force_over_length;
        auto& t_stiff = factory->get_stiffness();
        std::scoped_lock stiffness_lock(factory->get_stiffness_mutex());
        for(auto I = 0llu; I < chord.n_elem; ++I)
            for(auto J = 0llu; J < chord.n_elem; ++J) {
                t_stiff->at(dof_i(I), dof_i(J)) += der_repulsive(I, J);
                t_stiff->at(dof_j(I), dof_j(J)) += der_repulsive(I, J);
                t_stiff->at(dof_i(I), dof_j(J)) -= der_repulsive(I, J);
                t_stiff->at(dof_j(I), dof_i(J)) -= der_repulsive(I, J);
            }
    }
}

void HertzianDamped::apply(const shared_ptr<InteractionPair>& pair) const {
    const vec position_i = pair->position_i(), position_j = pair->position_j();
    const vec chord = position_i - position_j;
    const auto chord_length = norm(chord);
    const auto compression = pair->initial_gap() - chord_length;

    if(compression <= 0.) return;

    const vec unit_cord = chord / chord_length;
    const vec velocity_rel = pair->velocity_i() - pair->velocity_j();
    const auto velocity_projection = dot(velocity_rel, unit_cord);

    const auto pair_factor = four_third * std::sqrt(pair->effective_radius) * pair->effective_modulus * pair->effective_damping;
    const auto normal_force = pair_factor * std::pow(compression, .5) * velocity_projection;

    const uvec &dof_i = pair->dof_i(), &dof_j = pair->dof_j();

    {
        const vec repulsive = normal_force * unit_cord;
        auto& t_resistance = factory->modify_trial_constraint_resistance();
        std::scoped_lock resistance_lock(factory->get_trial_constraint_resistance_mutex());
        for(auto I = 0llu; I < chord.n_elem; ++I) {
            t_resistance(dof_i(I)) += repulsive(I);
            t_resistance(dof_j(I)) -= repulsive(I);
        }
    }

    {
        const mat der_unit_chord = (eye(chord.n_elem, chord.n_elem) - unit_cord * unit_cord.t()) / chord_length;
        const mat der_repulsive = normal_force * der_unit_chord + pair_factor * std::pow(compression, -.5) * unit_cord * (compression * velocity_rel.t() * der_unit_chord - .5 * velocity_projection * unit_cord.t());
        auto& t_stiff = factory->get_stiffness();
        std::scoped_lock stiffness_lock(factory->get_stiffness_mutex());
        for(auto I = 0llu; I < chord.n_elem; ++I)
            for(auto J = 0llu; J < chord.n_elem; ++J) {
                t_stiff->at(dof_i(I), dof_i(J)) += der_repulsive(I, J);
                t_stiff->at(dof_j(I), dof_j(J)) += der_repulsive(I, J);
                t_stiff->at(dof_i(I), dof_j(J)) -= der_repulsive(I, J);
                t_stiff->at(dof_j(I), dof_i(J)) -= der_repulsive(I, J);
            }
    }

    {
        const mat der_repulsive = pair_factor * std::pow(compression, .5) * unit_cord * unit_cord.t();
        auto& t_damping = factory->get_damping();
        std::scoped_lock damping_lock(factory->get_damping_mutex());
        for(auto I = 0llu; I < chord.n_elem; ++I)
            for(auto J = 0llu; J < chord.n_elem; ++J) {
                t_damping->at(dof_i(I), dof_i(J)) += der_repulsive(I, J);
                t_damping->at(dof_j(I), dof_j(J)) += der_repulsive(I, J);
                t_damping->at(dof_i(I), dof_j(J)) -= der_repulsive(I, J);
                t_damping->at(dof_j(I), dof_i(J)) -= der_repulsive(I, J);
            }
    }
}
