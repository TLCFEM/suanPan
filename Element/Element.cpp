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

#include "Element.h"
#include <Domain/DOF.h>
#include <Domain/DomainBase.h>
#include <Domain/Group/Group.h>
#include <Domain/Node.h>
#include <Material/Material.h>
#include <Section/Section.h>
#include <Toolbox/utility.h>

void Element::update_strain_energy() {
    if(trial_resistance.is_empty()) return;

    strain_energy += .5 * (current_resistance.is_empty() ? dot(get_incre_displacement(), trial_resistance) : dot(get_incre_displacement(), current_resistance + trial_resistance));
}

void Element::update_kinetic_energy() {
    if(trial_mass.is_empty()) return;

    const vec t_velocity = get_trial_velocity();

    kinetic_energy = .5 * dot(t_velocity, trial_mass * t_velocity);
}

void Element::update_viscous_energy() {
    if(trial_damping_force.is_empty()) return;

    viscous_energy += .5 * (current_damping_force.is_empty() ? dot(get_incre_displacement(), trial_damping_force) : dot(get_incre_displacement(), current_damping_force + trial_damping_force));
}

void Element::update_nonviscous_energy() {
    if(trial_nonviscous_force.is_empty()) return;

    nonviscous_energy += .5 * dot(get_incre_displacement(), real(sum(current_nonviscous_force + trial_nonviscous_force, 1)));
}

void Element::update_complementary_energy() {
    if(trial_resistance.is_empty()) return;

    const vec m_displacement = get_trial_displacement() + get_current_displacement();

    complementary_energy += .5 * (current_resistance.is_empty() ? dot(m_displacement, trial_resistance) : dot(m_displacement, trial_resistance - current_resistance));
}

void Element::update_momentum() {
    if(trial_mass.is_empty()) return;

    const vec t_velocity = get_trial_velocity();

    momentum = trial_mass * t_velocity;
}

/**
 * \brief generate a matrix that contains coordinates of connected nodes
 * \param num_dim number of dimension required
 * \return a matrix of following form
 *         | x_1  y_1  z_1  ... |
 *         | x_2  y_2  z_2  ... |
 *         | x_3  y_3  z_3  ... |
 */
mat Element::get_coordinate(const unsigned num_dim) const {
    mat ele_coor(num_node, num_dim, fill::zeros);

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_coor = node_ptr[I].lock()->get_coordinate();
        for(uword J = 0; J < std::min(static_cast<uword>(num_dim), t_coor.n_elem); ++J) ele_coor(I, J) = t_coor(J);
    }

    return ele_coor;
}

vec Element::get_incre_displacement() const {
    vec incre_displacement(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_disp = node_ptr[I].lock()->get_incre_displacement();
        for(unsigned J = 0; J < num_dof; ++J) incre_displacement(idx++) = t_disp(J);
    }

    return incre_displacement;
}

vec Element::get_incre_velocity() const {
    vec incre_velocity(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_vec = node_ptr[I].lock()->get_incre_velocity();
        for(unsigned J = 0; J < num_dof; ++J) incre_velocity(idx++) = t_vec(J);
    }

    return incre_velocity;
}

vec Element::get_incre_acceleration() const {
    vec incre_acceleration(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_acc = node_ptr[I].lock()->get_incre_acceleration();
        for(unsigned J = 0; J < num_dof; ++J) incre_acceleration(idx++) = t_acc(J);
    }

    return incre_acceleration;
}

vec Element::get_trial_displacement() const {
    vec trial_displacement(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_disp = node_ptr[I].lock()->get_trial_displacement();
        for(unsigned J = 0; J < num_dof; ++J) trial_displacement(idx++) = t_disp(J);
    }

    return trial_displacement;
}

vec Element::get_trial_velocity() const {
    vec trial_velocity(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_vec = node_ptr[I].lock()->get_trial_velocity();
        for(unsigned J = 0; J < num_dof; ++J) trial_velocity(idx++) = t_vec(J);
    }

    return trial_velocity;
}

vec Element::get_trial_acceleration() const {
    vec trial_acceleration(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_acc = node_ptr[I].lock()->get_trial_acceleration();
        for(unsigned J = 0; J < num_dof; ++J) trial_acceleration(idx++) = t_acc(J);
    }

    return trial_acceleration;
}

vec Element::get_current_displacement() const {
    vec current_displacement(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_disp = node_ptr[I].lock()->get_current_displacement();
        for(unsigned J = 0; J < num_dof; ++J) current_displacement(idx++) = t_disp(J);
    }

    return current_displacement;
}

vec Element::get_current_velocity() const {
    vec current_velocity(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_vec = node_ptr[I].lock()->get_current_velocity();
        for(unsigned J = 0; J < num_dof; ++J) current_velocity(idx++) = t_vec(J);
    }

    return current_velocity;
}

vec Element::get_current_acceleration() const {
    vec current_acceleration(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_acc = node_ptr[I].lock()->get_current_acceleration();
        for(unsigned J = 0; J < num_dof; ++J) current_acceleration(idx++) = t_acc(J);
    }

    return current_acceleration;
}

vec Element::get_node_incre_resistance() const {
    vec node_incre_resistance(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_force = node_ptr[I].lock()->get_incre_resistance();
        for(unsigned J = 0; J < num_dof; ++J) node_incre_resistance(idx++) = t_force(J);
    }

    return node_incre_resistance;
}

vec Element::get_node_trial_resistance() const {
    vec node_trial_resistance(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_force = node_ptr[I].lock()->get_trial_resistance();
        for(unsigned J = 0; J < num_dof; ++J) node_trial_resistance(idx++) = t_force(J);
    }

    return node_trial_resistance;
}

vec Element::get_node_current_resistance() const {
    vec node_current_resistance(num_size, fill::none);

    auto idx = 0;

    for(unsigned I = 0; I < num_node; ++I) {
        auto& t_force = node_ptr[I].lock()->get_current_resistance();
        for(unsigned J = 0; J < num_dof; ++J) node_current_resistance(idx++) = t_force(J);
    }

    return node_current_resistance;
}

std::vector<shared_ptr<Material>> Element::get_material(const shared_ptr<DomainBase>& D) const {
    std::vector<shared_ptr<Material>> material_pool;
    for(const auto& I : material_tag) material_pool.emplace_back(D->find<Material>(I) ? D->get<Material>(I) : nullptr);
    return material_pool;
}

std::vector<shared_ptr<Section>> Element::get_section(const shared_ptr<DomainBase>& D) const {
    std::vector<shared_ptr<Section>> section_pool;
    for(const auto& I : section_tag) section_pool.emplace_back(D->find<Section>(I) ? D->get<Section>(I) : nullptr);
    return section_pool;
}

Element::Element(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, std::vector<DOF>&& DI)
    : Element(T, NN, ND, std::forward<uvec>(NT), {}, false, MaterialType::D0, std::forward<std::vector<DOF>>(DI)) {}

Element::Element(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& MT, const bool F, const MaterialType MTP, std::vector<DOF>&& DI)
    : DataElement{std::forward<uvec>(NT), std::forward<uvec>(MT), uvec{}, F, true, true, true, true, true, {}}
    , ElementBase(T)
    , num_node(NN)
    , num_dof(ND)
    , mat_type(MTP)
    , sec_type(SectionType::D0)
    , dof_identifier(std::forward<std::vector<DOF>>(DI)) { suanpan_assert([&] { if(!dof_identifier.empty() && num_dof != dof_identifier.size()) throw invalid_argument("size of dof identifier must meet number of dofs"); }); }

Element::Element(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& ST, const bool F, const SectionType STP, std::vector<DOF>&& DI)
    : DataElement{std::forward<uvec>(NT), uvec{}, std::forward<uvec>(ST), F, true, true, true, true, true, {}}
    , ElementBase(T)
    , num_node(NN)
    , num_dof(ND)
    , mat_type(MaterialType::D0)
    , sec_type(STP)
    , dof_identifier(std::forward<std::vector<DOF>>(DI)) { suanpan_assert([&] { if(!dof_identifier.empty() && num_dof != dof_identifier.size()) throw invalid_argument("size of dof identifier must meet number of dofs"); }); }

// for contact elements that use node groups
Element::Element(const unsigned T, const unsigned ND, uvec&& GT)
    : DataElement{std::forward<uvec>(GT), {}, {}, false, true, true, true, true, true, {}}
    , ElementBase(T)
    , num_node(static_cast<unsigned>(-1))
    , num_dof(ND)
    , use_group(true)
    , mat_type(MaterialType::D0)
    , sec_type(SectionType::D0) {}

// for elements that use other elements
Element::Element(const unsigned T, const unsigned ND, const unsigned ET, const unsigned NT)
    : DataElement{{NT}, {}, {}, false, true, true, true, true, true, {}}
    , ElementBase(T)
    , num_node(static_cast<unsigned>(-1))
    , num_dof(ND)
    , use_other(ET)
    , mat_type(MaterialType::D0)
    , sec_type(SectionType::D0) {}

int Element::initialize_base(const shared_ptr<DomainBase>& D) {
    // initialized already, check node validity
    if(node_ptr.size() == num_node) {
        for(const auto& I : node_ptr) if(const auto t_node = I.lock(); nullptr == t_node || !t_node->is_active()) return SUANPAN_FAIL;
        return SUANPAN_SUCCESS;
    }

    // use node group instead of node
    if(use_group) {
        std::vector<const uvec*> pool;
        pool.reserve(node_encoding.n_elem);
        for(const auto I : node_encoding)
            if(D->find<Group>(I)) pool.emplace_back(&D->get<Group>(I)->get_pool());
            else return SUANPAN_FAIL;

        uword counter = 0;
        for(const auto& I : pool) counter += I->size();
        access::rw(num_node) = static_cast<unsigned>(counter);

        auto& n_encoding = access::rw(node_encoding);
        n_encoding.zeros(num_node);
        counter = 0;
        for(const auto& I : pool) for(const auto J : *I) n_encoding(counter++) = J;
    }

    // embedded elements use other elements
    if(0 != use_other) {
        if(!D->find<Element>(use_other)) return SUANPAN_FAIL;

        unsigned size = 1;
        auto& t_element = D->get<Element>(use_other);
        size += t_element->get_node_number();
        access::rw(num_node) = size;

        auto& n_encoding = access::rw(node_encoding);
        n_encoding.resize(size);
        n_encoding.tail(size - 1llu) = t_element->get_node_encoding();
    }

    // first initialization
    access::rw(num_size) = num_node * num_dof;

    if(0 == num_size) return SUANPAN_FAIL;

    dof_encoding.set_size(num_size);

    // check if nodes are still valid
    auto inactive = false;
    node_ptr.clear();
    node_ptr.reserve(num_node);
    for(const auto& t_tag : node_encoding) {
        if(!D->find<Node>(t_tag)) {
            suanpan_warning("Element {} disabled as node {} cannot be found.\n", get_tag(), t_tag);
            return SUANPAN_FAIL;
        }
        auto& t_node = D->get<Node>(t_tag);
        node_ptr.emplace_back(t_node);
        if(!t_node->is_active()) inactive = true;
    }
    if(inactive) {
        suanpan_warning("Element {} disabled as inactive nodes used.\n", get_tag());
        return SUANPAN_FAIL;
    }
    for(auto& I : node_ptr) if(I.lock()->get_dof_number() < num_dof) I.lock()->set_dof_number(num_dof);

    // check if material models are valid
    if(MaterialType::D0 != mat_type)
        for(const auto& t_tag : material_tag)
            if(auto& t_material = D->get<Material>(t_tag); nullptr == t_material || !t_material->is_active() || t_material->get_material_type() != MaterialType::DS && t_material->get_material_type() != mat_type) {
                suanpan_warning("Element {} disabled as material {} cannot be found or type mismatch.\n", get_tag(), t_tag);
                return SUANPAN_FAIL;
            }

    // check if section models are valid
    if(SectionType::D0 != sec_type)
        for(const auto& t_tag : section_tag)
            if(auto& t_section = D->get<Section>(t_tag); nullptr == t_section || !t_section->is_active() || t_section->get_section_type() != sec_type) {
                suanpan_warning("Element {} disabled as section {} cannot be found or type mismatch.\n", get_tag(), t_tag);
                return SUANPAN_FAIL;
            }

#ifdef SUANPAN_VTK
    // vtk visualization setup
    Setup();
#endif

    return SUANPAN_SUCCESS;
}

void Element::set_initialized(const bool F) const { access::rw(initialized) = F; }

void Element::set_symmetric(const bool F) const { access::rw(symmetric) = F; }

bool Element::is_initialized() const { return initialized; }

bool Element::is_symmetric() const { return symmetric; }

bool Element::is_nlgeom() const { return nlgeom; }

void Element::update_dof_encoding() {
    unsigned idx = 0;
    for(const auto& tmp_ptr : node_ptr) {
        auto& node_dof = tmp_ptr.lock()->get_reordered_dof();
        for(unsigned i = 0; i < num_dof; ++i) dof_encoding(idx++) = node_dof(i);
    }

    dof_mapping.clear();
    dof_mapping.reserve(num_size);
    const uvec dof_index = sort_index(dof_encoding), dof_reordered = dof_encoding(dof_index);
    for(auto I = 0llu; I < dof_index.n_elem; ++I) for(auto J = I; J < dof_index.n_elem; ++J) dof_mapping.emplace_back(MappingDOF{dof_reordered(J), dof_reordered(I), dof_index(J), dof_index(I)});
    std::sort(dof_mapping.begin(), dof_mapping.end(), [](const MappingDOF& A, const MappingDOF& B) { return A.l_col == B.l_col ? A.l_row < B.l_row : A.l_col < B.l_col; });

    if(!dof_identifier.empty()) for(const auto& tmp_ptr : node_ptr) tmp_ptr.lock()->set_dof_identifier(dof_identifier);
}

bool Element::if_update_mass() const { return update_mass; }

bool Element::if_update_damping() const { return update_damping; }

bool Element::if_update_nonviscous() const { return update_nonviscous; }

bool Element::if_update_stiffness() const { return update_stiffness; }

bool Element::if_update_geometry() const { return update_geometry; }

const uvec& Element::get_dof_encoding() const { return dof_encoding; }

const uvec& Element::get_node_encoding() const { return node_encoding; }

const std::vector<MappingDOF>& Element::get_dof_mapping() const { return dof_mapping; }

const uvec& Element::get_material_tag() const { return material_tag; }

const uvec& Element::get_section_tag() const { return section_tag; }

unsigned Element::get_dof_number() const { return num_dof; }

unsigned Element::get_node_number() const { return num_node; }

unsigned Element::get_total_number() const { return num_size; }

void Element::clear_node_ptr() { node_ptr.clear(); }

const std::vector<weak_ptr<Node>>& Element::get_node_ptr() const { return node_ptr; }

const vec& Element::get_trial_resistance() const { return trial_resistance; }

const vec& Element::get_current_resistance() const { return current_resistance; }

const vec& Element::get_trial_damping_force() const { return trial_damping_force; }

const vec& Element::get_current_damping_force() const { return current_damping_force; }

const cx_mat& Element::get_trial_nonviscous_force() const { return trial_nonviscous_force; }

const cx_mat& Element::get_current_nonviscous_force() const { return current_nonviscous_force; }

const vec& Element::get_trial_inertial_force() {
    if(!trial_mass.empty()) trial_inertial_force = trial_mass * get_trial_acceleration();
    return trial_inertial_force;
}

const vec& Element::get_current_inertial_force() {
    if(!current_mass.empty()) current_inertial_force = current_mass * get_current_acceleration();
    return current_inertial_force;
}

const vec& Element::get_trial_body_force() const { return trial_body_force; }

const vec& Element::get_current_body_force() const { return current_body_force; }

const vec& Element::get_trial_traction() const { return trial_traction; }

const vec& Element::get_current_traction() const { return current_traction; }

const mat& Element::get_trial_mass() const { return trial_mass; }

const mat& Element::get_trial_damping() const { return trial_damping; }

const mat& Element::get_trial_nonviscous() const { return trial_nonviscous; }

const mat& Element::get_trial_stiffness() const { return trial_stiffness; }

const mat& Element::get_trial_geometry() const { return trial_geometry; }

const mat& Element::get_trial_secant() const { return get_trial_stiffness(); }

const mat& Element::get_current_mass() const { return current_mass; }

const mat& Element::get_current_damping() const { return current_damping; }

const mat& Element::get_current_nonviscous() const { return current_nonviscous; }

const mat& Element::get_current_stiffness() const { return current_stiffness; }

const mat& Element::get_current_geometry() const { return current_geometry; }

const mat& Element::get_current_secant() const { return get_current_stiffness(); }

const mat& Element::get_initial_mass() const { return initial_mass; }

const mat& Element::get_initial_damping() const { return initial_damping; }

const mat& Element::get_initial_nonviscous() const { return initial_nonviscous; }

const mat& Element::get_initial_stiffness() const { return initial_stiffness; }

const mat& Element::get_initial_geometry() const { return initial_geometry; }

const mat& Element::get_initial_secant() const { return get_initial_stiffness(); }

const mat& Element::get_mass_container() const { return mass_container; }

const mat& Element::get_stiffness_container() const { return stiffness_container; }

int Element::clear_status() {
    if(update_mass) trial_mass = current_mass = initial_mass;
    if(update_damping) trial_damping = current_damping = initial_damping;
    if(update_nonviscous) trial_nonviscous = current_nonviscous = initial_nonviscous;
    if(update_stiffness) trial_stiffness = current_stiffness = initial_stiffness;
    if(update_geometry) trial_geometry = current_geometry = initial_geometry;

    if(!trial_resistance.is_empty()) trial_resistance.zeros();
    if(!current_resistance.is_empty()) current_resistance.zeros();
    if(!trial_damping_force.is_empty()) trial_damping_force.zeros();
    if(!current_damping_force.is_empty()) current_damping_force.zeros();
    if(!trial_nonviscous_force.is_empty()) trial_nonviscous_force.zeros();
    if(!current_nonviscous_force.is_empty()) current_nonviscous_force.zeros();
    if(!trial_inertial_force.is_empty()) trial_inertial_force.zeros();
    if(!current_inertial_force.is_empty()) current_inertial_force.zeros();

    // do not clear node pointer here, it will be cleared in global initialization
    // node_ptr.clear();

    strain_energy = 0.;
    kinetic_energy = 0.;
    viscous_energy = 0.;
    nonviscous_energy = 0.;
    complementary_energy = 0.;
    momentum.zeros();

    set_initialized(false);

    return SUANPAN_SUCCESS;
}

int Element::commit_status() {
    update_strain_energy();
    update_kinetic_energy();
    update_viscous_energy();
    update_nonviscous_energy();
    update_complementary_energy();
    update_momentum();

    if(update_mass && !trial_mass.is_empty()) current_mass = trial_mass;
    if(update_damping && !trial_damping.is_empty()) current_damping = trial_damping;
    if(update_stiffness && !trial_stiffness.is_empty()) current_stiffness = trial_stiffness;
    if(update_geometry && !trial_geometry.is_empty()) current_geometry = trial_geometry;
    if(!trial_resistance.is_empty()) current_resistance = trial_resistance;
    if(!trial_damping_force.is_empty()) current_damping_force = trial_damping_force;
    if(!trial_nonviscous_force.is_empty()) current_nonviscous_force = trial_nonviscous_force;
    if(!trial_inertial_force.is_empty()) current_inertial_force = trial_inertial_force;

    return SUANPAN_SUCCESS;
}

int Element::reset_status() {
    if(update_mass && !trial_mass.is_empty()) trial_mass = current_mass;
    if(update_damping && !trial_damping.is_empty()) trial_damping = current_damping;
    if(update_nonviscous && !trial_nonviscous.is_empty()) trial_nonviscous = current_nonviscous;
    if(update_stiffness && !trial_stiffness.is_empty()) trial_stiffness = current_stiffness;
    if(update_geometry && !trial_geometry.is_empty()) trial_geometry = current_geometry;
    if(!trial_resistance.is_empty()) trial_resistance = current_resistance;
    if(!trial_damping_force.is_empty()) trial_damping_force = current_damping_force;
    if(!trial_nonviscous_force.is_empty()) trial_nonviscous_force = current_nonviscous_force;
    if(!trial_inertial_force.is_empty()) trial_inertial_force = current_inertial_force;

    return SUANPAN_SUCCESS;
}

const vec& Element::update_body_force(const vec& load_factor) { return body_force.is_empty() ? trial_body_force : trial_body_force = body_force * load_factor; }

const vec& Element::update_traction(const vec& load_factor) { return traction.is_empty() ? trial_traction : trial_traction = traction * load_factor; }

std::vector<vec> Element::record(const OutputType) { return {}; }

double Element::get_strain_energy() const { return strain_energy; }

double Element::get_complementary_energy() const { return complementary_energy; }

double Element::get_kinetic_energy() const { return kinetic_energy; }

double Element::get_viscous_energy() const { return viscous_energy; }

double Element::get_nonviscous_energy() const { return nonviscous_energy; }

const vec& Element::get_momentum() const { return momentum; }

double Element::get_momentum_component(const DOF D) const {
    auto [flag, position] = if_contain(dof_identifier, D);

    if(!flag || momentum.empty()) return 0.;

    auto momentum_component = 0.;
    for(auto I = 0u; I < num_node; ++I, position += num_dof) momentum_component += momentum(position);

    return momentum_component;
}

double Element::get_characteristic_length() const { return characteristic_length; }

mat Element::compute_shape_function(const mat&, unsigned) const { return {}; }

void ConstantMass(DataElement* E) {
    E->update_mass = false;
    E->current_mass = mat(E->initial_mass.memptr(), E->initial_mass.n_rows, E->initial_mass.n_cols, false, true);
    E->trial_mass = mat(E->initial_mass.memptr(), E->initial_mass.n_rows, E->initial_mass.n_cols, false, true);
}

void ConstantDamping(DataElement* E) {
    E->update_damping = false;
    E->current_damping = mat(E->initial_damping.memptr(), E->initial_damping.n_rows, E->initial_damping.n_cols, false, true);
    E->trial_damping = mat(E->initial_damping.memptr(), E->initial_damping.n_rows, E->initial_damping.n_cols, false, true);
}

void ConstantStiffness(DataElement* E) {
    E->update_stiffness = false;
    E->current_stiffness = mat(E->initial_stiffness.memptr(), E->initial_stiffness.n_rows, E->initial_stiffness.n_cols, false, true);
    E->trial_stiffness = mat(E->initial_stiffness.memptr(), E->initial_stiffness.n_rows, E->initial_stiffness.n_cols, false, true);
}

void ConstantGeometry(DataElement* E) {
    E->update_geometry = false;
    E->current_geometry = mat(E->initial_geometry.memptr(), E->initial_geometry.n_rows, E->initial_geometry.n_cols, false, true);
    E->trial_geometry = mat(E->initial_geometry.memptr(), E->initial_geometry.n_rows, E->initial_geometry.n_cols, false, true);
}

mat get_coordinate(const ElementBase* const E, const unsigned N) { return E->get_coordinate(N); }

vec get_incre_displacement(const ElementBase* const E) { return E->get_incre_displacement(); }

vec get_incre_velocity(const ElementBase* const E) { return E->get_incre_velocity(); }

vec get_incre_acceleration(const ElementBase* const E) { return E->get_incre_acceleration(); }

vec get_trial_displacement(const ElementBase* const E) { return E->get_trial_displacement(); }

vec get_trial_velocity(const ElementBase* const E) { return E->get_trial_velocity(); }

vec get_trial_acceleration(const ElementBase* const E) { return E->get_trial_acceleration(); }

vec get_current_displacement(const ElementBase* const E) { return E->get_current_displacement(); }

vec get_current_velocity(const ElementBase* const E) { return E->get_current_velocity(); }

vec get_current_acceleration(const ElementBase* const E) { return E->get_current_acceleration(); }
