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

#include "Domain.h"
#include <Constraint/Constraint.h>
#include <Constraint/Criterion/Criterion.h>
#include <Converger/Converger.h>
#include <Database/Database.h>
#include <Domain/FactoryHelper.hpp>
#include <Domain/Group.h>
#include <Domain/Node.h>
#include <Element/Element.h>
#include <Element/Modifier/Modifier.h>
#include <Element/Utility/Orientation.h>
#include <Load/Amplitude/Amplitude.h>
#include <Load/Load.h>
#include <Material/Material.h>
#include <Recorder/Recorder.h>
#include <Section/Section.h>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Solver.h>
#include <Step/Step.h>
#include <Toolbox/sort_color.hpp>
#include <Toolbox/sort_rcm.h>
#include <numeric>

#ifdef SUANPAN_HDF5
#include <hdf5/hdf5.h>
#include <hdf5/hdf5_hl.h>
#endif

Domain::Domain(const unsigned T)
    : DomainBase(T)
    , factory(make_shared<Factory<double>>()) {}

Domain::~Domain() { for(const auto& I : thread_pond) I->get(); }

void Domain::set_factory(const shared_ptr<LongFactory>& F) {
    if(factory == F) return;

    factory = F;
    updated = false;
}

const shared_ptr<LongFactory>& Domain::get_factory() const { return factory; }

bool Domain::insert(const shared_ptr<future<void>>& T) {
    thread_pond.emplace_back(T);
    return true;
}

void Domain::wait() { for(const auto& thread : thread_pond) thread->wait(); }

bool Domain::insert(const shared_ptr<ExternalModule>& E) {
    external_module_pond.emplace_back(E);
    return true;
}

const ExternalModuleQueue& Domain::get_external_module_pool() const { return external_module_pond; }

bool Domain::insert(const shared_ptr<Amplitude>& A) {
    updated = false;
    return amplitude_pond.insert(A);
}

bool Domain::insert(const shared_ptr<Constraint>& C) {
    updated = false;
    return constraint_pond.insert(C);
}

bool Domain::insert(const shared_ptr<Converger>& C) {
    updated = false;
    return converger_pond.insert(C);
}

bool Domain::insert(const shared_ptr<Criterion>& C) {
    updated = false;
    return criterion_pond.insert(C);
}

bool Domain::insert(const shared_ptr<Database>& D) {
    updated = false;
    return database_pond.insert(D);
}

bool Domain::insert(const shared_ptr<Element>& E) {
    updated = false;
    return element_pond.insert(E);
}

bool Domain::insert(const shared_ptr<Group>& E) {
    updated = false;
    return group_pond.insert(E);
}

bool Domain::insert(const shared_ptr<Integrator>& I) {
    updated = false;
    return integrator_pond.insert(I);
}

bool Domain::insert(const shared_ptr<Load>& L) {
    updated = false;
    return load_pond.insert(L);
}

bool Domain::insert(const shared_ptr<Material>& M) {
    updated = false;
    return material_pond.insert(M);
}

bool Domain::insert(const shared_ptr<Modifier>& M) {
    updated = false;
    return modifier_pond.insert(M);
}

bool Domain::insert(const shared_ptr<Node>& N) {
    updated = false;
    return node_pond.insert(N);
}

bool Domain::insert(const shared_ptr<Orientation>& N) {
    updated = false;
    return orientation_pond.insert(N);
}

bool Domain::insert(const shared_ptr<Recorder>& R) {
    updated = false;
    return recorder_pond.insert(R);
}

bool Domain::insert(const shared_ptr<Section>& S) {
    updated = false;
    return section_pond.insert(S);
}

bool Domain::insert(const shared_ptr<Solver>& S) {
    updated = false;
    return solver_pond.insert(S);
}

bool Domain::insert(const shared_ptr<Step>& S) {
    updated = false;
    return step_pond.insert({S->get_tag(), S}).second;
}

bool Domain::erase_amplitude(const unsigned T) {
    if(!find<Amplitude>(T)) return true;

    updated = false;
    return amplitude_pond.erase(T);
}

bool Domain::erase_constraint(const unsigned T) {
    if(!find<Constraint>(T)) return true;

    updated = false;
    return constraint_pond.erase(T);
}

bool Domain::erase_converger(const unsigned T) {
    if(!find<Converger>(T)) return true;

    updated = false;
    return converger_pond.erase(T);
}

bool Domain::erase_criterion(const unsigned T) {
    if(!find<Criterion>(T)) return true;

    updated = false;
    return criterion_pond.erase(T);
}

bool Domain::erase_database(const unsigned T) {
    if(!find<Database>(T)) return true;

    updated = false;
    return database_pond.erase(T);
}

bool Domain::erase_element(const unsigned T) {
    if(!find<Element>(T)) return true;

    updated = false;
    return element_pond.erase(T);
}

bool Domain::erase_group(const unsigned T) {
    if(!find<Group>(T)) return true;

    updated = false;
    return group_pond.erase(T);
}

bool Domain::erase_integrator(const unsigned T) {
    if(!find<Integrator>(T)) return true;

    updated = false;
    return integrator_pond.erase(T);
}

bool Domain::erase_load(const unsigned T) {
    if(!find<Load>(T)) return true;

    updated = false;
    return load_pond.erase(T);
}

bool Domain::erase_material(const unsigned T) {
    if(!find<Material>(T)) return true;

    updated = false;
    return material_pond.erase(T);
}

bool Domain::erase_modifier(const unsigned T) {
    if(!find<Modifier>(T)) return true;

    updated = false;
    return modifier_pond.erase(T);
}

bool Domain::erase_node(const unsigned T) {
    if(!find<Node>(T)) return true;

    updated = false;
    return node_pond.erase(T);
}

bool Domain::erase_orientation(const unsigned T) {
    if(!find<Orientation>(T)) return true;

    updated = false;
    return orientation_pond.erase(T);
}

bool Domain::erase_recorder(const unsigned T) {
    if(!find<Recorder>(T)) return true;

    updated = false;
    return recorder_pond.erase(T);
}

bool Domain::erase_section(const unsigned T) {
    if(!find<Section>(T)) return true;

    updated = false;
    return section_pond.erase(T);
}

bool Domain::erase_solver(const unsigned T) {
    if(!find<Solver>(T)) return true;

    updated = false;
    return solver_pond.erase(T);
}

bool Domain::erase_step(const unsigned T) {
    if(!find<Step>(T)) return true;

    updated = false;
    return step_pond.erase(T) == 1;
}

void Domain::disable_amplitude(const unsigned T) {
    if(!find<Amplitude>(T) || !get<Amplitude>(T)->is_active() || get<Amplitude>(T)->is_guarded()) return;

    updated = false;
    amplitude_pond.disable(T);
}

void Domain::disable_constraint(const unsigned T) {
    if(!find<Constraint>(T) || !get<Constraint>(T)->is_active() || get<Constraint>(T)->is_guarded()) return;

    updated = false;
    constraint_pond.disable(T);
}

void Domain::disable_converger(const unsigned T) {
    if(!find<Converger>(T) || !get<Converger>(T)->is_active() || get<Converger>(T)->is_guarded()) return;

    updated = false;
    converger_pond.disable(T);
}

void Domain::disable_criterion(const unsigned T) {
    if(!find<Criterion>(T) || !get<Criterion>(T)->is_active() || get<Criterion>(T)->is_guarded()) return;

    updated = false;
    criterion_pond.disable(T);
}

void Domain::disable_database(const unsigned T) {
    if(!find<Database>(T) || !get<Database>(T)->is_active() || get<Database>(T)->is_guarded()) return;

    updated = false;
    database_pond.disable(T);
}

void Domain::disable_element(const unsigned T) {
    if(!find<Element>(T) || !get<Element>(T)->is_active() || get<Element>(T)->is_guarded()) return;

    updated = false;
    element_pond.disable(T);
}

void Domain::disable_group(const unsigned T) {
    if(!find<Group>(T) || !get<Group>(T)->is_active() || get<Group>(T)->is_guarded()) return;

    updated = false;
    group_pond.disable(T);
}

void Domain::disable_integrator(const unsigned T) {
    if(!find<Integrator>(T) || !get<Integrator>(T)->is_active() || get<Integrator>(T)->is_guarded()) return;

    updated = false;
    integrator_pond.disable(T);
}

void Domain::disable_load(const unsigned T) {
    if(!find<Load>(T) || !get<Load>(T)->is_active() || get<Load>(T)->is_guarded()) return;

    updated = false;
    load_pond.disable(T);
}

void Domain::disable_material(const unsigned T) {
    if(!find<Material>(T) || !get<Material>(T)->is_active() || get<Material>(T)->is_guarded()) return;

    updated = false;
    material_pond.disable(T);
}

void Domain::disable_modifier(const unsigned T) {
    if(!find<Modifier>(T) || !get<Modifier>(T)->is_active() || get<Modifier>(T)->is_guarded()) return;

    updated = false;
    modifier_pond.disable(T);
}

void Domain::disable_node(const unsigned T) {
    if(!find<Node>(T) || !get<Node>(T)->is_active() || get<Node>(T)->is_guarded()) return;

    updated = false;
    node_pond.disable(T);
}

void Domain::disable_orientation(const unsigned T) {
    if(!find<Orientation>(T) || !get<Orientation>(T)->is_active() || get<Orientation>(T)->is_guarded()) return;

    updated = false;
    orientation_pond.disable(T);
}

void Domain::disable_recorder(const unsigned T) {
    if(!find<Recorder>(T) || !get<Recorder>(T)->is_active() || get<Recorder>(T)->is_guarded()) return;

    updated = false;
    recorder_pond.disable(T);
}

void Domain::disable_section(const unsigned T) {
    if(!find<Section>(T) || !get<Section>(T)->is_active() || get<Section>(T)->is_guarded()) return;

    updated = false;
    section_pond.disable(T);
}

void Domain::disable_solver(const unsigned T) {
    if(!find<Solver>(T) || !get<Solver>(T)->is_active() || get<Solver>(T)->is_guarded()) return;

    updated = false;
    solver_pond.disable(T);
}

void Domain::disable_step(const unsigned T) {
    if(!find<Step>(T) || !get<Step>(T)->is_active() || get<Step>(T)->is_guarded()) return;

    updated = false;

    step_pond.at(T)->disable();
}

void Domain::enable_amplitude(const unsigned T) {
    if(!find<Amplitude>(T) || get<Amplitude>(T)->is_active()) return;

    updated = false;
    amplitude_pond.enable(T);
}

void Domain::enable_constraint(const unsigned T) {
    if(!find<Constraint>(T) || get<Constraint>(T)->is_active()) return;

    updated = false;
    constraint_pond.enable(T);
}

void Domain::enable_converger(const unsigned T) {
    if(!find<Converger>(T) || get<Converger>(T)->is_active()) return;

    updated = false;
    converger_pond.enable(T);
}

void Domain::enable_criterion(const unsigned T) {
    if(!find<Criterion>(T) || get<Criterion>(T)->is_active()) return;

    updated = false;
    criterion_pond.enable(T);
}

void Domain::enable_database(const unsigned T) {
    if(!find<Database>(T) || get<Database>(T)->is_active()) return;

    updated = false;
    database_pond.enable(T);
}

void Domain::enable_element(const unsigned T) {
    if(!find<Element>(T) || get<Element>(T)->is_active()) return;

    updated = false;
    element_pond.enable(T);
}

void Domain::enable_group(const unsigned T) {
    if(!find<Group>(T) || get<Group>(T)->is_active()) return;

    updated = false;
    group_pond.enable(T);
}

void Domain::enable_integrator(const unsigned T) {
    if(!find<Integrator>(T) || get<Integrator>(T)->is_active()) return;

    updated = false;
    integrator_pond.enable(T);
}

void Domain::enable_load(const unsigned T) {
    if(!find<Load>(T) || get<Load>(T)->is_active()) return;

    updated = false;
    load_pond.enable(T);
}

void Domain::enable_material(const unsigned T) {
    if(!find<Material>(T) || get<Material>(T)->is_active()) return;

    updated = false;
    material_pond.enable(T);
}

void Domain::enable_modifier(const unsigned T) {
    if(!find<Modifier>(T) || get<Modifier>(T)->is_active()) return;

    updated = false;
    modifier_pond.enable(T);
}

void Domain::enable_node(const unsigned T) {
    if(!find<Node>(T) || get<Node>(T)->is_active()) return;

    updated = false;
    node_pond.enable(T);
}

void Domain::enable_orientation(const unsigned T) {
    if(!find<Orientation>(T) || get<Orientation>(T)->is_active()) return;

    updated = false;
    orientation_pond.enable(T);
}

void Domain::enable_recorder(const unsigned T) {
    if(!find<Recorder>(T) || get<Recorder>(T)->is_active()) return;

    updated = false;
    recorder_pond.enable(T);
}

void Domain::enable_section(const unsigned T) {
    if(!find<Section>(T) || get<Section>(T)->is_active()) return;

    updated = false;
    section_pond.enable(T);
}

void Domain::enable_solver(const unsigned T) {
    if(!find<Solver>(T) || get<Solver>(T)->is_active()) return;

    updated = false;
    solver_pond.enable(T);
}

void Domain::enable_step(const unsigned T) {
    if(!find<Step>(T) || get<Step>(T)->is_active()) return;

    updated = false;

    step_pond.at(T)->enable();
}

const shared_ptr<Amplitude>& Domain::get_amplitude(const unsigned T) const { return amplitude_pond.at(T); }

const shared_ptr<Constraint>& Domain::get_constraint(const unsigned T) const { return constraint_pond.at(T); }

const shared_ptr<Converger>& Domain::get_converger(const unsigned T) const { return converger_pond.at(T); }

const shared_ptr<Criterion>& Domain::get_criterion(const unsigned T) const { return criterion_pond.at(T); }

const shared_ptr<Database>& Domain::get_database(const unsigned T) const { return database_pond.at(T); }

const shared_ptr<Element>& Domain::get_element(const unsigned T) const { return element_pond.at(T); }

const shared_ptr<Group>& Domain::get_group(const unsigned T) const { return group_pond.at(T); }

const shared_ptr<Integrator>& Domain::get_integrator(const unsigned T) const { return integrator_pond.at(T); }

const shared_ptr<Load>& Domain::get_load(const unsigned T) const { return load_pond.at(T); }

const shared_ptr<Material>& Domain::get_material(const unsigned T) const { return material_pond.at(T); }

const shared_ptr<Modifier>& Domain::get_modifier(const unsigned T) const { return modifier_pond.at(T); }

const shared_ptr<Node>& Domain::get_node(const unsigned T) const { return node_pond.at(T); }

const shared_ptr<Orientation>& Domain::get_orientation(const unsigned T) const { return orientation_pond.at(T); }

const shared_ptr<Recorder>& Domain::get_recorder(const unsigned T) const { return recorder_pond.at(T); }

const shared_ptr<Section>& Domain::get_section(const unsigned T) const { return section_pond.at(T); }

const shared_ptr<Solver>& Domain::get_solver(const unsigned T) const { return solver_pond.at(T); }

const shared_ptr<Step>& Domain::get_step(const unsigned T) const { return step_pond.at(T); }

const AmplitudeQueue& Domain::get_amplitude_pool() const { return amplitude_pond.get(); }

const ConstraintQueue& Domain::get_constraint_pool() const { return constraint_pond.get(); }

const ConvergerQueue& Domain::get_converger_pool() const { return converger_pond.get(); }

const CriterionQueue& Domain::get_criterion_pool() const { return criterion_pond.get(); }

const DatabaseQueue& Domain::get_database_pool() const { return database_pond.get(); }

const ElementQueue& Domain::get_element_pool() const { return element_pond.get(); }

const GroupQueue& Domain::get_group_pool() const { return group_pond.get(); }

const IntegratorQueue& Domain::get_integrator_pool() const { return integrator_pond.get(); }

const LoadQueue& Domain::get_load_pool() const { return load_pond.get(); }

const MaterialQueue& Domain::get_material_pool() const { return material_pond.get(); }

const ModifierQueue& Domain::get_modifier_pool() const { return modifier_pond.get(); }

const NodeQueue& Domain::get_node_pool() const { return node_pond.get(); }

const OrientationQueue& Domain::get_orientation_pool() const { return orientation_pond.get(); }

const RecorderQueue& Domain::get_recorder_pool() const { return recorder_pond.get(); }

const SectionQueue& Domain::get_section_pool() const { return section_pond.get(); }

const SolverQueue& Domain::get_solver_pool() const { return solver_pond.get(); }

const StepQueue& Domain::get_step_pool() const { return step_pond; }

size_t Domain::get_amplitude() const { return amplitude_pond.size(); }

size_t Domain::get_constraint() const { return constraint_pond.size(); }

size_t Domain::get_converger() const { return converger_pond.size(); }

size_t Domain::get_criterion() const { return criterion_pond.size(); }

size_t Domain::get_database() const { return database_pond.size(); }

size_t Domain::get_element() const { return element_pond.size(); }

size_t Domain::get_group() const { return group_pond.size(); }

size_t Domain::get_integrator() const { return integrator_pond.size(); }

size_t Domain::get_load() const { return load_pond.size(); }

size_t Domain::get_material() const { return material_pond.size(); }

size_t Domain::get_modifier() const { return modifier_pond.size(); }

size_t Domain::get_node() const { return node_pond.size(); }

size_t Domain::get_orientation() const { return orientation_pond.size(); }

size_t Domain::get_recorder() const { return recorder_pond.size(); }

size_t Domain::get_section() const { return section_pond.size(); }

size_t Domain::get_solver() const { return solver_pond.size(); }

size_t Domain::get_step() const { return step_pond.size(); }

bool Domain::find_amplitude(const unsigned T) const { return amplitude_pond.find(T); }

bool Domain::find_constraint(const unsigned T) const { return constraint_pond.find(T); }

bool Domain::find_converger(const unsigned T) const { return converger_pond.find(T); }

bool Domain::find_criterion(const unsigned T) const { return criterion_pond.find(T); }

bool Domain::find_database(const unsigned T) const { return database_pond.find(T); }

bool Domain::find_element(const unsigned T) const { return element_pond.find(T); }

bool Domain::find_group(const unsigned T) const { return group_pond.find(T); }

bool Domain::find_integrator(const unsigned T) const { return integrator_pond.find(T); }

bool Domain::find_load(const unsigned T) const { return load_pond.find(T); }

bool Domain::find_material(const unsigned T) const { return material_pond.find(T); }

bool Domain::find_modifier(const unsigned T) const { return modifier_pond.find(T); }

bool Domain::find_node(const unsigned T) const { return node_pond.find(T); }

bool Domain::find_orientation(const unsigned T) const { return orientation_pond.find(T); }

bool Domain::find_recorder(const unsigned T) const { return recorder_pond.find(T); }

bool Domain::find_section(const unsigned T) const { return section_pond.find(T); }

bool Domain::find_solver(const unsigned T) const { return solver_pond.find(T); }

bool Domain::find_step(const unsigned T) const { return step_pond.contains(T); }

void Domain::set_current_step_tag(const unsigned T) { current_step_tag = T; }

void Domain::set_current_converger_tag(const unsigned T) { current_converger_tag = T; }

void Domain::set_current_integrator_tag(const unsigned T) { current_integrator_tag = T; }

void Domain::set_current_solver_tag(const unsigned T) { current_solver_tag = T; }

unsigned Domain::get_current_step_tag() { return current_step_tag; }

unsigned Domain::get_current_converger_tag() { return current_converger_tag; }

unsigned Domain::get_current_integrator_tag() { return current_integrator_tag; }

unsigned Domain::get_current_solver_tag() { return current_solver_tag; }

const shared_ptr<Step>& Domain::get_current_step() const { return get_step(current_step_tag); }

const shared_ptr<Converger>& Domain::get_current_converger() const { return get_converger(current_converger_tag); }

const shared_ptr<Integrator>& Domain::get_current_integrator() const { return get_integrator(current_integrator_tag); }

const shared_ptr<Solver>& Domain::get_current_solver() const { return get_solver(current_solver_tag); }

void Domain::insert_loaded_dof(const uvec& T) { loaded_dofs.insert(T.cbegin(), T.cend()); }

void Domain::insert_restrained_dof(const uvec& T) { restrained_dofs.insert(T.cbegin(), T.cend()); }

void Domain::insert_constrained_dof(const uvec& T) { constrained_dofs.insert(T.cbegin(), T.cend()); }

void Domain::insert_loaded_dof(const uword T) { loaded_dofs.insert(T); }

void Domain::insert_restrained_dof(const uword T) { restrained_dofs.insert(T); }

void Domain::insert_constrained_dof(const uword T) { constrained_dofs.insert(T); }

const suanpan::unordered_set<uword>& Domain::get_loaded_dof() const { return loaded_dofs; }

const suanpan::unordered_set<uword>& Domain::get_restrained_dof() const { return restrained_dofs; }

const suanpan::unordered_set<uword>& Domain::get_constrained_dof() const { return constrained_dofs; }

bool Domain::is_updated() const { return updated.load(); }

bool Domain::is_sparse() const { return factory->is_sparse(); }

void Domain::set_color_model(const ColorMethod B) {
    color_model = B;
    color_map.clear();
}

const vector<vector<unsigned>>& Domain::get_color_map() const { return color_map; }

/**
 * \brief list of connected element tags for each node
 * \return list
 */
vector<vector<uword>> Domain::get_node_connectivity() {
    unsigned max_tag = 0;
    for(const auto& [t_tag , t_node] : node_pond) if(t_tag > max_tag) max_tag = t_tag;

    vector node_connectivity(++max_tag, vector<uword>{});

    for(const auto& [t_tag, t_element] : element_pond) for(const auto t_node : t_element->get_node_encoding()) node_connectivity[t_node].emplace_back(t_tag);

    suanpan_for_each(node_connectivity.begin(), node_connectivity.end(), [](vector<uword>& t_node) { suanpan::unique(t_node); });

    return node_connectivity;
}

/**
 * \brief list of connected node tags for each element
 * \return list
 */
vector<uvec> Domain::get_element_connectivity() {
    unsigned max_tag = 0;
    for(const auto& [t_tag, t_element] : element_pond) if(t_tag > max_tag) max_tag = t_tag;

    vector element_connectivity(++max_tag, uvec{});

    suanpan_for_each(element_pond.cbegin(), element_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Element>>& t_element) { element_connectivity[t_element.first] = t_element.second->get_node_encoding(); });

    return element_connectivity;
}

int Domain::reorder_dof() {
    // assign dof label for active dof
    unsigned dof_counter = 0;
    for(const auto& t_node : node_pond.get()) t_node->set_original_dof(dof_counter);
    if(0 == dof_counter) return SUANPAN_FAIL;
    // active flag is now properly set for node and element

    // RCM optimization
    // collect connectivity
#ifdef SUANPAN_MT
    vector<concurrent_unordered_set<uword>> adjacency(dof_counter);
#else
    vector<unordered_set<uword>> adjacency(dof_counter);
#endif
    suanpan_for_each(element_pond.get().cbegin(), element_pond.get().cend(), [&](const shared_ptr<Element>& t_element) {
        t_element->update_dof_encoding();
        for(auto& t_encoding = t_element->get_dof_encoding(); const auto& I : t_encoding) for(const auto& J : t_encoding) adjacency[I].insert(J);
    });

    // for nonlinear constraint
    for(const auto& [t_tag, t_constraint] : constraint_pond) {
        if(!t_constraint->is_connected()) continue;
        vector<uword> t_encoding;
        for(auto& I : t_constraint->get_node_encoding()) if(find<Node>(I)) for(auto& J : get<Node>(I)->get_reordered_dof()) t_encoding.emplace_back(J);
        for(const auto& I : t_encoding) for(const auto& J : t_encoding) adjacency[I].insert(J);
    }

    // count number of degree
    uvec num_degree(dof_counter);
    suanpan_for(0u, dof_counter, [&](const unsigned I) { num_degree(I) = adjacency[I].size(); });

    // sort each column according to its degree
    vector<uvec> adjacency_sorted(dof_counter);
    suanpan_for(0u, dof_counter, [&](const unsigned i) {
        uvec t_vec(num_degree(i));
        unsigned j = 0;
        for(const auto k : adjacency[i]) t_vec(j++) = k;
        adjacency_sorted[i] = t_vec(sort_index(num_degree(t_vec)));
    });

    const auto idx_rcm = sort_rcm(adjacency_sorted, num_degree);
    const uvec idx_sorted = sort_index(idx_rcm);

    // get bandwidth
    auto low_bw = 0, up_bw = 0;
    for(unsigned i = 0; i < dof_counter; ++i)
        for(const auto j : adjacency[idx_rcm(i)]) {
            if(const auto t_bw = static_cast<int>(idx_sorted(j)) - static_cast<int>(i); t_bw > low_bw) low_bw = t_bw;
            else if(t_bw < up_bw) up_bw = t_bw;
        }

    // assign new labels to active nodes
    auto& t_node_pond = node_pond.get();
    suanpan_for_each(t_node_pond.cbegin(), t_node_pond.cend(), [&](const shared_ptr<Node>& t_node) { t_node->set_reordered_dof(idx_sorted(t_node->get_original_dof())); });

    factory->set_size(dof_counter);
    factory->set_bandwidth(static_cast<unsigned>(low_bw), static_cast<unsigned>(-up_bw));

    return SUANPAN_SUCCESS;
}

int Domain::assign_color() {
    // deal with k-coloring optimization
    if(ColorMethod::OFF != color_model) {
        const auto color_algorithm = ColorMethod::WP == color_model ? sort_color_wp<unsigned> : sort_color_mis<unsigned>;

        auto node_tag = 0llu;
        std::unordered_map<uword, uword> node_map; // old_tag -> new_tag
        for(auto& t_node : node_pond.get()) node_map[t_node->get_tag()] = node_tag++;

        auto element_tag = 0u;
        vector<unsigned> element_map; // new_idx -> old_tag
        element_map.reserve(element_pond.get().size());
        vector<vector<unsigned>> node_register(node_tag);
        for(auto& t_element : element_pond.get()) {
            element_map.emplace_back(t_element->get_tag());
            for(const auto t_node : t_element->get_node_encoding()) node_register[node_map.at(t_node)].emplace_back(element_tag);
            element_tag++;
        }

        suanpan::graph<unsigned> element_register(element_tag);

        suanpan_for_each(node_register.begin(), node_register.end(), [&](const vector<unsigned>& node) { for(const auto I : node) for(const auto J : node) element_register[I].insert(J); });

        color_map = color_algorithm(element_register);

        suanpan_for_each(color_map.begin(), color_map.end(), [&](vector<unsigned>& color) { std::ranges::transform(color, color.begin(), [&](const unsigned element) { return element_map[element]; }); });

        suanpan_debug("The model is colored by %llu colors.\n", color_map.size());
    }

    // count how many entries in the sparse form and preallocate memory
    if(is_sparse()) {
        auto& t_element_pool = element_pond.get();
        factory->set_entry(std::transform_reduce(t_element_pool.cbegin(), t_element_pool.cend(), 1000llu, std::plus(), [](const shared_ptr<Element>& t_element) { return t_element->get_total_number() * t_element->get_total_number(); }));
    }

    return SUANPAN_SUCCESS;
}

int Domain::restart() {
    // try to initialize to check if anything changes
    if(SUANPAN_SUCCESS != initialize()) return SUANPAN_FAIL;

    // try to initialize to check if anything changes
    if(SUANPAN_SUCCESS != factory->initialize()) return SUANPAN_FAIL;

    // initialize load and constraint for current step
    if(SUANPAN_SUCCESS != initialize_load()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != initialize_constraint()) return SUANPAN_FAIL;

    // check-in existing nodal quantities
    if(SUANPAN_SUCCESS != update_current_status()) return SUANPAN_FAIL;

    // reference may depend on reformulated displacement
    if(SUANPAN_SUCCESS != initialize_reference()) return SUANPAN_FAIL;

    return SUANPAN_SUCCESS;
}

int Domain::soft_restart() { return updated ? SUANPAN_SUCCESS : restart(); }

int Domain::initialize() {
    if(updated) return SUANPAN_SUCCESS;

    // some independent objects
    converger_pond.update();
    integrator_pond.update();
    solver_pond.update();

    // for restart analysis
    suanpan_for_each(load_pond.cbegin(), load_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Load>>& t_load) { t_load.second->set_initialized(false); });
    suanpan_for_each(constraint_pond.cbegin(), constraint_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Constraint>>& t_constraint) { t_constraint.second->set_initialized(false); });

    // amplitude should be updated before load
    suanpan_for_each(amplitude_pond.cbegin(), amplitude_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Amplitude>>& t_amplitude) { t_amplitude.second->initialize(shared_from_this()); });
    amplitude_pond.update();

    suanpan::set<unsigned> remove_list;

    remove_list.clear();
    suanpan_for_each(material_pond.cbegin(), material_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Material>>& t_material) {
        if(t_material.second->is_initialized() || !t_material.second->is_active()) return;

        // if fail to initialize, the material is invalid, remove it from the model
        if(SUANPAN_SUCCESS != t_material.second->initialize_base(shared_from_this()) || SUANPAN_SUCCESS != t_material.second->initialize(shared_from_this())) {
            disable_material(t_material.first);
            remove_list.insert(t_material.first);
            return;
        }

        t_material.second->set_initialized(true);
    });
    for(const auto I : remove_list) erase_material(I);
    material_pond.update();

    remove_list.clear();
    suanpan_for_each(section_pond.cbegin(), section_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Section>>& t_section) {
        if(t_section.second->is_initialized() || !t_section.second->is_active()) return;

        // if fail to initialize, the section is invalid, remove it from the model
        if(SUANPAN_SUCCESS != t_section.second->initialize_base(shared_from_this()) || SUANPAN_SUCCESS != t_section.second->initialize(shared_from_this())) {
            disable_section(t_section.first);
            remove_list.insert(t_section.first);
            return;
        }

        t_section.second->set_initialized(true);
    });
    for(const auto I : remove_list) erase_section(I);
    section_pond.update();

    // set dof number to zero before first initialisation of elements
    suanpan_for_each(node_pond.cbegin(), node_pond.cend(), [](const std::pair<unsigned, shared_ptr<Node>>& t_node) {
        // for restart analysis that may reassign dof
        t_node.second->set_initialized(false);
        t_node.second->set_dof_number(0);
    });
    node_pond.update();

    // groups reply on nodes
    suanpan_for_each(group_pond.cbegin(), group_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Group>>& t_group) { t_group.second->initialize(shared_from_this()); });
    group_pond.update();

    // element may reply on groups
    remove_list.clear();
    suanpan_for_each(element_pond.cbegin(), element_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Element>>& t_element) {
        if(!t_element.second->is_active()) return;
        // clear node pointer array to enable initialisation
        t_element.second->clear_node_ptr();
        if(SUANPAN_SUCCESS != t_element.second->initialize_base(shared_from_this())) {
            disable_element(t_element.first);
            remove_list.insert(t_element.first);
        }
    });
    for(const auto I : remove_list) erase_element(I);
    element_pond.update();

    // initialise nodal variables with proper number of DoFs
    suanpan_for_each(node_pond.cbegin(), node_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Node>>& t_node) { t_node.second->initialize(shared_from_this()); });
    node_pond.update();

    // order matters between element base initialization and derived method call
    updated = true;

    if(SUANPAN_SUCCESS != reorder_dof()) return SUANPAN_FAIL;

    remove_list.clear();
    auto nlgeom = false;
    // initialize derived elements
    auto& t_element_pond = element_pond.get();
    suanpan_for_each(t_element_pond.cbegin(), t_element_pond.cend(), [&](const shared_ptr<Element>& t_element) {
        if(t_element->is_initialized()) {
            // model changed due to addition/deletion of elements, etc.
            // the element itself is not modified, just update dof encoding
            t_element->update_dof_encoding();
            // cannot update status here due to some unknown reason
            // maybe machine error or inconsistent tolerance due to material models
            if(t_element->is_nlgeom()) nlgeom = true;
            return;
        }

        // if first initialisation fails, the element can be safely deleted
        if(SUANPAN_SUCCESS != t_element->initialize(shared_from_this())) {
            disable_element(t_element->get_tag());
            remove_list.insert(t_element->get_tag());
            return;
        }

        t_element->set_initialized(true);
        t_element->update_dof_encoding();
        // to get the very first matrix
        t_element->update_status();
        if(t_element->is_nlgeom()) nlgeom = true;
    });
    for(const auto I : remove_list) erase_element(I);
    element_pond.update();
    if(t_element_pond.empty()) {
        suanpan_warning("no active element.\n");
        return SUANPAN_FAIL;
    }

    // set some factory properties for later initialisation
    factory->set_nlgeom(nlgeom);

    // initialize modifier based on updated element pool
    suanpan_for_each(modifier_pond.cbegin(), modifier_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Modifier>>& t_modifier) { t_modifier.second->initialize(shared_from_this()); });
    modifier_pond.update();
    // sort to ensure lower performs first
    if(auto& t_modifier_pool = access::rw(modifier_pond.get()); t_modifier_pool.size() > 1)
        suanpan_sort(t_modifier_pool.begin(), t_modifier_pool.end(), [&](const shared_ptr<Modifier>& a, const shared_ptr<Modifier>& b) { return a->get_tag() < b->get_tag(); });

    // recorder may depend on groups, nodes, elements, etc.
    suanpan_for_each(recorder_pond.cbegin(), recorder_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Recorder>>& t_recorder) { t_recorder.second->initialize(shared_from_this()); });
    recorder_pond.update();

    suanpan_for_each(criterion_pond.cbegin(), criterion_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Criterion>>& t_criterion) { t_criterion.second->initialize(shared_from_this()); });
    criterion_pond.update();

    // element initialization may change the status of the domain
    return updated ? assign_color() : initialize();
}

int Domain::initialize_load() {
    suanpan_for_each(load_pond.cbegin(), load_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Load>>& t_load) { if(t_load.second->validate_step(shared_from_this()) && !t_load.second->is_initialized() && SUANPAN_FAIL == t_load.second->initialize(shared_from_this())) disable_load(t_load.first); });
    load_pond.update();

    factory->update_reference_size();

    return SUANPAN_SUCCESS;
}

int Domain::initialize_constraint() {
    suanpan_for_each(constraint_pond.cbegin(), constraint_pond.cend(), [&](const std::pair<unsigned, shared_ptr<Constraint>>& t_constraint) { if(t_constraint.second->validate_step(shared_from_this()) && !t_constraint.second->is_initialized() && SUANPAN_FAIL == t_constraint.second->initialize(shared_from_this())) disable_constraint(t_constraint.first); });
    constraint_pond.update();

    return SUANPAN_SUCCESS;
}

int Domain::initialize_reference() {
    // used in displacement controlled algorithm
    // need to initialize it for updating of load factor
    const auto ref_dof = to_uvec(factory->get_reference_dof());

    if(ref_dof.is_empty()) return SUANPAN_SUCCESS;

    factory->initialize_settlement();
    // check in current displacement in case it is nonzero
    auto ref_settlement = factory->get_trial_settlement();
    ref_settlement(ref_dof) = factory->get_trial_displacement()(ref_dof);
    factory->set_trial_settlement(ref_settlement);
    factory->set_current_settlement(ref_settlement);
    factory->initialize_load_factor();
    auto& ref_load = get_reference_load(factory);
    for(uword I = 0; I < ref_dof.n_elem; ++I) ref_load(ref_dof(I), I) = 1.;

    return SUANPAN_SUCCESS;
}

int Domain::process_load(const bool full) {
    loaded_dofs.clear();

    auto& t_load = get_trial_load(factory);
    if(!t_load.empty()) t_load.zeros();

    auto& t_settlement = get_trial_settlement(factory);
    if(!t_settlement.empty()) t_settlement.zeros();

    const auto process_handler = full ? std::mem_fn(&Load::process) : std::mem_fn(&Load::process_resistance);

    auto code = 0;
    for(auto& I : load_pond.get()) {
        if(!I->is_initialized()) continue;
        code += std::invoke(process_handler, I, shared_from_this());
        if(!I->get_trial_load().empty()) t_load += I->get_trial_load();
        if(!I->get_trial_settlement().empty()) t_settlement += I->get_trial_settlement();
    }

    factory->update_trial_load(t_load);
    factory->update_trial_settlement(t_settlement);

    return code;
}

int Domain::process_constraint(const bool full) {
    // ! Tasks:
    // ! 1. process initialised constraints
    // ! 2. clear previous auxiliary variables as they may change size
    // ! 3. assemble constraint resistance
    // ! 4. assemble auxiliary resistance, load, and stiffness
    // ! 5. record auxiliary encoding
    // ! assemble stiffness shall be done in integrator

    restrained_dofs.clear();
    constrained_dofs.clear();

    auto& constraint_resistance = get_trial_constraint_resistance(factory);
    constraint_resistance.zeros();

    const auto process_handler = full ? std::mem_fn(&Constraint::process) : std::mem_fn(&Constraint::process_resistance);

    auto code = 0;
    auto counter = 0u;
    for(auto& I : constraint_pond.get()) {
        if(!I->is_initialized()) continue;
        code += std::invoke(process_handler, I, shared_from_this());
        if(!I->get_resistance().empty()) constraint_resistance += I->get_resistance();
        counter += I->get_multiplier_size();
    }

    factory->set_mpc(counter);

    if(0 == counter) return code;

    auto& t_encoding = get_auxiliary_encoding(factory);
    auto& t_resistance = get_auxiliary_resistance(factory);
    auto& t_load = get_auxiliary_load(factory);
    auto& t_stiffness = get_auxiliary_stiffness(factory);

    counter = 0u;
    for(auto& I : constraint_pond.get()) {
        const auto m_size = I->get_multiplier_size();
        if(!I->is_initialized() || 0 == m_size) continue;
        const auto e_size = counter + m_size - 1;
        t_encoding.subvec(counter, e_size).fill(I->get_tag());
        if(!I->get_auxiliary_resistance().empty()) t_resistance.subvec(counter, e_size) = I->get_auxiliary_resistance();
        if(!I->get_auxiliary_load().empty()) t_load.subvec(counter, e_size) = I->get_auxiliary_load();
        if(!I->get_auxiliary_stiffness().empty()) t_stiffness.cols(counter, e_size) = I->get_auxiliary_stiffness();
        counter += m_size;
    }

    return code;
}

int Domain::process_criterion() {
    auto code = 0;
    for(auto& I : criterion_pond.get()) code += I->process(shared_from_this());
    return code;
}

int Domain::process_modifier() {
    auto code = 0;
    // use sequential for_each only
    for(auto& I : modifier_pond.get()) code += I->update_status();
    return code;
}

void Domain::record() {
    auto& t_recorder_pool = recorder_pond.get();
    suanpan_for_each(t_recorder_pool.cbegin(), t_recorder_pool.cend(), [&](const shared_ptr<Recorder>& t_recorder) { t_recorder->record(shared_from_this()); });
}

void Domain::enable_all() {
    updated = false;

    amplitude_pond.enable();
    constraint_pond.enable();
    converger_pond.enable();
    criterion_pond.enable();
    element_pond.enable();
    integrator_pond.enable();
    load_pond.enable();
    material_pond.enable();
    node_pond.enable();
    recorder_pond.enable();
    section_pond.enable();
    solver_pond.enable();
}

void Domain::summary() const {
    suanpan_info("Domain %u contains:\n\t%llu nodes, %llu elements, %llu materials,\n", get_tag(), get_node(), get_element(), get_material());
    suanpan_info("\t%llu loads, %llu constraints and %llu recorders.\n", get_load(), get_constraint(), get_recorder());
}

void Domain::erase_machine_error() const {
    auto& t_ninja = get_ninja(factory);
    for(const auto& I : restrained_dofs) t_ninja(I) = 0.;
}

void Domain::update_load() {}

void Domain::update_constraint() {
    auto& t_encoding = factory->get_auxiliary_encoding();
    auto& t_lambda = factory->get_auxiliary_lambda();

    for(auto I = 0llu; I < t_encoding.n_elem;) {
        auto& t_constraint = get<Constraint>(t_encoding(I));
        const auto t_size = t_constraint->get_multiplier_size();
        t_constraint->update_status(t_lambda.subvec(I, I + t_size - 1));
        I += t_size;
    }
}

void Domain::assemble_load_stiffness() {}

void Domain::assemble_constraint_stiffness() { for(auto& I : get_constraint_pool()) if(I->is_initialized() && !I->get_stiffness().empty()) factory->assemble_stiffness(I->get_stiffness(), I->get_dof_encoding()); }

void Domain::save(string file_name) {
#ifdef SUANPAN_HDF5
    if(string::npos == file_name.find(".h5") && string::npos == file_name.find(".H5")) file_name += ".h5";

    const auto file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    const string group_name = "/Node";

    const auto group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for(const auto& t_node : get_node_pool()) {
        const auto data = t_node->save();

        const hsize_t dimension[2] = {data.n_rows, data.n_cols};

        H5LTmake_dataset(group_id, std::to_string(t_node->get_tag()).c_str(), 2, dimension, H5T_NATIVE_DOUBLE, data.mem);
    }

    H5Gclose(group_id);
    H5Fclose(file_id);
#endif
}
