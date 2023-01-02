/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "Rebar2D.h"
#include <Domain/DomainBase.h>
#include <Toolbox/tensorToolbox.h>

Rebar2D::Rebar2D(const unsigned T, const unsigned XT, const unsigned YT, const double RX, const double RY)
    : Material2D(T, PlaneType::S, 0.)
    , tag_x(XT)
    , tag_y(YT)
    , ratio_x(RX)
    , ratio_y(RY) {}

Rebar2D::Rebar2D(const Rebar2D& old_obj)
    : Material2D(old_obj)
    , tag_x(old_obj.tag_x)
    , tag_y(old_obj.tag_y)
    , ratio_x(old_obj.ratio_x)
    , ratio_y(old_obj.ratio_y)
    , rebar_x(suanpan::make_copy(old_obj.rebar_x))
    , rebar_y(suanpan::make_copy(old_obj.rebar_y)) {}

int Rebar2D::initialize(const shared_ptr<DomainBase>& D) {
    rebar_x = suanpan::initialized_material_copy(D, tag_x);
    rebar_y = suanpan::initialized_material_copy(D, tag_y);

    if(nullptr == rebar_x || nullptr == rebar_y || rebar_x->get_material_type() != MaterialType::D1 || rebar_y->get_material_type() != MaterialType::D1) {
        SP_E("A valid 1D host material is required.\n");
        return SUANPAN_FAIL;
    }

    access::rw(density) = ratio_x * rebar_x->get_parameter(ParameterType::DENSITY) + ratio_y * rebar_y->get_parameter(ParameterType::DENSITY);

    initial_stiffness.zeros(3, 3);

    initial_stiffness(0, 0) = ratio_x * rebar_x->get_initial_stiffness().at(0);
    initial_stiffness(1, 1) = ratio_y * rebar_y->get_initial_stiffness().at(0);

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Rebar2D::get_copy() { return make_unique<Rebar2D>(*this); }

int Rebar2D::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;

    // update status
    if(SUANPAN_SUCCESS != rebar_x->update_trial_status(trial_strain(0))) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != rebar_y->update_trial_status(trial_strain(1))) return SUANPAN_FAIL;

    trial_stress.set_size(3);

    // collect main stress components
    trial_stress(0) = ratio_x * rebar_x->get_trial_stress().at(0);
    trial_stress(1) = ratio_y * rebar_y->get_trial_stress().at(0);
    trial_stress(2) = 0.;

    // collect principal stiffness components
    trial_stiffness.zeros(3, 3);
    trial_stiffness(0, 0) = ratio_x * rebar_x->get_trial_stiffness().at(0);
    trial_stiffness(1, 1) = ratio_y * rebar_y->get_trial_stiffness().at(0);

    return SUANPAN_SUCCESS;
}

int Rebar2D::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    return rebar_x->clear_status() + rebar_y->clear_status();
}

int Rebar2D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return rebar_x->commit_status() + rebar_y->commit_status();
}

int Rebar2D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return rebar_x->reset_status() + rebar_y->reset_status();
}

vector<vec> Rebar2D::record(const OutputType P) {
    vector<vec> data;

    for(const auto& I : rebar_x->record(P)) data.emplace_back(I);
    for(const auto& I : rebar_y->record(P)) data.emplace_back(I);

    return data;
}

void Rebar2D::print() {
    suanpan_info("A rebar layer with major/minor reinforcement ratio of %.3E and %.3E.\n", ratio_x, ratio_y);
    suanpan_info("Major: ");
    if(rebar_x) rebar_x->print();
    suanpan_info("Minor: ");
    if(rebar_y) rebar_y->print();
}
