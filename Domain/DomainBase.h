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
/**
 * @class DomainBase
 * @brief The DomainBase class is a template.
 *
 * The DomainBase is simply an abstraction of the Domain class. It provides all
 * methods signature that are used in Domain class. The purpose is to split the
 * declaration and implementation apart. As the Domain class is widely used in
 * many other classes. The dependency hierarchy is simplified if replaced by the
 * DomainBase.
 *
 * @author tlc
 * @date 01/10/2017
 * @version 0.2.0
 * @file DomainBase.h
 * @addtogroup Domain
 * @{
 */

#ifndef DOMAINBASE_H
#define DOMAINBASE_H

#include <future>
#include <Domain/Tag.h>
#include <Toolbox/container.h>

using std::future;

template<sp_d T> class Factory;
class Amplitude;
class Constraint;
class Converger;
class Criterion;
class Database;
class Element;
class ExternalModule;
class Group;
class Integrator;
class Load;
class Material;
class Modifier;
class Node;
class Orientation;
class Recorder;
class Section;
class Solver;
class Step;

using AmplitudeQueue = std::vector<shared_ptr<Amplitude>>;
using ConstraintQueue = std::vector<shared_ptr<Constraint>>;
using ConvergerQueue = std::vector<shared_ptr<Converger>>;
using CriterionQueue = std::vector<shared_ptr<Criterion>>;
using DatabaseQueue = std::vector<shared_ptr<Database>>;
using ElementQueue = std::vector<shared_ptr<Element>>;
using GroupQueue = std::vector<shared_ptr<Group>>;
using IntegratorQueue = std::vector<shared_ptr<Integrator>>;
using LoadQueue = std::vector<shared_ptr<Load>>;
using MaterialQueue = std::vector<shared_ptr<Material>>;
using ModifierQueue = std::vector<shared_ptr<Modifier>>;
using NodeQueue = std::vector<shared_ptr<Node>>;
using OrientationQueue = std::vector<shared_ptr<Orientation>>;
using RecorderQueue = std::vector<shared_ptr<Recorder>>;
using SectionQueue = std::vector<shared_ptr<Section>>;
using SolverQueue = std::vector<shared_ptr<Solver>>;
using StepQueue = std::map<unsigned, shared_ptr<Step>>;

using LongFactory = Factory<double>;

enum class ColorMethod {
    OFF,
    WP,
    MIS
};

class DomainBase : public Tag {
public:
    explicit DomainBase(unsigned);
    DomainBase(const DomainBase&) = delete;            // copy forbidden
    DomainBase(DomainBase&&) = delete;                 // move forbidden
    DomainBase& operator=(const DomainBase&) = delete; // assign forbidden
    DomainBase& operator=(DomainBase&&) = delete;      // assign forbidden
    ~DomainBase() override;

    virtual void set_factory(const shared_ptr<LongFactory>&) = 0;
    [[nodiscard]] virtual const shared_ptr<LongFactory>& get_factory() const = 0;

    virtual bool insert(const shared_ptr<future<void>>&) = 0;

    virtual void wait() = 0;

    virtual bool insert(const shared_ptr<ExternalModule>&) = 0;
    [[nodiscard]] virtual const std::vector<shared_ptr<ExternalModule>>& get_external_module_pool() const = 0;

    virtual bool insert(const shared_ptr<Amplitude>&) = 0;
    virtual bool insert(const shared_ptr<Constraint>&) = 0;
    virtual bool insert(const shared_ptr<Converger>&) = 0;
    virtual bool insert(const shared_ptr<Criterion>&) = 0;
    virtual bool insert(const shared_ptr<Database>&) = 0;
    virtual bool insert(const shared_ptr<Element>&) = 0;
    virtual bool insert(const shared_ptr<Group>&) = 0;
    virtual bool insert(const shared_ptr<Integrator>&) = 0;
    virtual bool insert(const shared_ptr<Load>&) = 0;
    virtual bool insert(const shared_ptr<Material>&) = 0;
    virtual bool insert(const shared_ptr<Modifier>&) = 0;
    virtual bool insert(const shared_ptr<Node>&) = 0;
    virtual bool insert(const shared_ptr<Orientation>&) = 0;
    virtual bool insert(const shared_ptr<Recorder>&) = 0;
    virtual bool insert(const shared_ptr<Section>&) = 0;
    virtual bool insert(const shared_ptr<Solver>&) = 0;
    virtual bool insert(const shared_ptr<Step>&) = 0;

    template<typename T> bool erase(unsigned);
    virtual bool erase_amplitude(unsigned) = 0;
    virtual bool erase_constraint(unsigned) = 0;
    virtual bool erase_converger(unsigned) = 0;
    virtual bool erase_criterion(unsigned) = 0;
    virtual bool erase_database(unsigned) = 0;
    virtual bool erase_element(unsigned) = 0;
    virtual bool erase_group(unsigned) = 0;
    virtual bool erase_integrator(unsigned) = 0;
    virtual bool erase_load(unsigned) = 0;
    virtual bool erase_material(unsigned) = 0;
    virtual bool erase_modifier(unsigned) = 0;
    virtual bool erase_node(unsigned) = 0;
    virtual bool erase_orientation(unsigned) = 0;
    virtual bool erase_recorder(unsigned) = 0;
    virtual bool erase_section(unsigned) = 0;
    virtual bool erase_solver(unsigned) = 0;
    virtual bool erase_step(unsigned) = 0;

    virtual void disable_amplitude(unsigned) = 0;
    virtual void disable_constraint(unsigned) = 0;
    virtual void disable_converger(unsigned) = 0;
    virtual void disable_criterion(unsigned) = 0;
    virtual void disable_database(unsigned) = 0;
    virtual void disable_element(unsigned) = 0;
    virtual void disable_group(unsigned) = 0;
    virtual void disable_integrator(unsigned) = 0;
    virtual void disable_load(unsigned) = 0;
    virtual void disable_material(unsigned) = 0;
    virtual void disable_modifier(unsigned) = 0;
    virtual void disable_node(unsigned) = 0;
    virtual void disable_orientation(unsigned) = 0;
    virtual void disable_recorder(unsigned) = 0;
    virtual void disable_section(unsigned) = 0;
    virtual void disable_solver(unsigned) = 0;
    virtual void disable_step(unsigned) = 0;

    virtual void enable_amplitude(unsigned) = 0;
    virtual void enable_constraint(unsigned) = 0;
    virtual void enable_converger(unsigned) = 0;
    virtual void enable_criterion(unsigned) = 0;
    virtual void enable_database(unsigned) = 0;
    virtual void enable_element(unsigned) = 0;
    virtual void enable_group(unsigned) = 0;
    virtual void enable_integrator(unsigned) = 0;
    virtual void enable_load(unsigned) = 0;
    virtual void enable_material(unsigned) = 0;
    virtual void enable_modifier(unsigned) = 0;
    virtual void enable_node(unsigned) = 0;
    virtual void enable_orientation(unsigned) = 0;
    virtual void enable_recorder(unsigned) = 0;
    virtual void enable_section(unsigned) = 0;
    virtual void enable_solver(unsigned) = 0;
    virtual void enable_step(unsigned) = 0;

    template<typename T> const shared_ptr<T>& get(unsigned);
    template<typename T> const shared_ptr<T>& get(uword);
    template<typename T> std::vector<shared_ptr<T>> get(const uvec&);
    [[nodiscard]] virtual const shared_ptr<Amplitude>& get_amplitude(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Constraint>& get_constraint(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Converger>& get_converger(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Criterion>& get_criterion(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Database>& get_database(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Element>& get_element(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Group>& get_group(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Integrator>& get_integrator(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Load>& get_load(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Material>& get_material(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Modifier>& get_modifier(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Node>& get_node(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Orientation>& get_orientation(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Recorder>& get_recorder(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Section>& get_section(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Solver>& get_solver(unsigned) const = 0;
    [[nodiscard]] virtual const shared_ptr<Step>& get_step(unsigned) const = 0;

    template<typename T> const std::vector<shared_ptr<T>>& get_pool();
    [[nodiscard]] virtual const AmplitudeQueue& get_amplitude_pool() const = 0;
    [[nodiscard]] virtual const ConstraintQueue& get_constraint_pool() const = 0;
    [[nodiscard]] virtual const ConvergerQueue& get_converger_pool() const = 0;
    [[nodiscard]] virtual const CriterionQueue& get_criterion_pool() const = 0;
    [[nodiscard]] virtual const DatabaseQueue& get_database_pool() const = 0;
    [[nodiscard]] virtual const ElementQueue& get_element_pool() const = 0;
    [[nodiscard]] virtual const GroupQueue& get_group_pool() const = 0;
    [[nodiscard]] virtual const IntegratorQueue& get_integrator_pool() const = 0;
    [[nodiscard]] virtual const LoadQueue& get_load_pool() const = 0;
    [[nodiscard]] virtual const MaterialQueue& get_material_pool() const = 0;
    [[nodiscard]] virtual const ModifierQueue& get_modifier_pool() const = 0;
    [[nodiscard]] virtual const NodeQueue& get_node_pool() const = 0;
    [[nodiscard]] virtual const OrientationQueue& get_orientation_pool() const = 0;
    [[nodiscard]] virtual const RecorderQueue& get_recorder_pool() const = 0;
    [[nodiscard]] virtual const SectionQueue& get_section_pool() const = 0;
    [[nodiscard]] virtual const SolverQueue& get_solver_pool() const = 0;
    [[nodiscard]] virtual const StepQueue& get_step_pool() const = 0;

    friend shared_ptr<Amplitude>& get_amplitude(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Constraint>& get_constraint(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Converger>& get_converger(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Criterion>& get_criterion(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Database>& get_database(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Element>& get_element(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Group>& get_group(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Integrator>& get_integrator(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Load>& get_load(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Material>& get_material(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Modifier>& get_modifier(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Node>& get_node(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Orientation>& get_orientation(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Recorder>& get_recorder(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Section>& get_section(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Solver>& get_solver(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Step>& get_step(const shared_ptr<DomainBase>&, unsigned);

    template<typename T> size_t get();
    [[nodiscard]] virtual size_t get_amplitude() const = 0;
    [[nodiscard]] virtual size_t get_constraint() const = 0;
    [[nodiscard]] virtual size_t get_converger() const = 0;
    [[nodiscard]] virtual size_t get_criterion() const = 0;
    [[nodiscard]] virtual size_t get_database() const = 0;
    [[nodiscard]] virtual size_t get_element() const = 0;
    [[nodiscard]] virtual size_t get_group() const = 0;
    [[nodiscard]] virtual size_t get_integrator() const = 0;
    [[nodiscard]] virtual size_t get_load() const = 0;
    [[nodiscard]] virtual size_t get_material() const = 0;
    [[nodiscard]] virtual size_t get_modifier() const = 0;
    [[nodiscard]] virtual size_t get_node() const = 0;
    [[nodiscard]] virtual size_t get_orientation() const = 0;
    [[nodiscard]] virtual size_t get_recorder() const = 0;
    [[nodiscard]] virtual size_t get_section() const = 0;
    [[nodiscard]] virtual size_t get_solver() const = 0;
    [[nodiscard]] virtual size_t get_step() const = 0;

    template<typename T> bool find(unsigned);
    template<typename T> bool find(uword);
    template<typename T> bool find(const uvec&);
    [[nodiscard]] virtual bool find_amplitude(unsigned) const = 0;
    [[nodiscard]] virtual bool find_constraint(unsigned) const = 0;
    [[nodiscard]] virtual bool find_converger(unsigned) const = 0;
    [[nodiscard]] virtual bool find_criterion(unsigned) const = 0;
    [[nodiscard]] virtual bool find_database(unsigned) const = 0;
    [[nodiscard]] virtual bool find_element(unsigned) const = 0;
    [[nodiscard]] virtual bool find_group(unsigned) const = 0;
    [[nodiscard]] virtual bool find_integrator(unsigned) const = 0;
    [[nodiscard]] virtual bool find_load(unsigned) const = 0;
    [[nodiscard]] virtual bool find_material(unsigned) const = 0;
    [[nodiscard]] virtual bool find_modifier(unsigned) const = 0;
    [[nodiscard]] virtual bool find_node(unsigned) const = 0;
    [[nodiscard]] virtual bool find_orientation(unsigned) const = 0;
    [[nodiscard]] virtual bool find_recorder(unsigned) const = 0;
    [[nodiscard]] virtual bool find_section(unsigned) const = 0;
    [[nodiscard]] virtual bool find_solver(unsigned) const = 0;
    [[nodiscard]] virtual bool find_step(unsigned) const = 0;

    virtual void set_current_step_tag(unsigned) = 0;
    virtual void set_current_converger_tag(unsigned) = 0;
    virtual void set_current_integrator_tag(unsigned) = 0;
    virtual void set_current_solver_tag(unsigned) = 0;

    virtual unsigned get_current_step_tag() = 0;
    virtual unsigned get_current_converger_tag() = 0;
    virtual unsigned get_current_integrator_tag() = 0;
    virtual unsigned get_current_solver_tag() = 0;

    [[nodiscard]] virtual const shared_ptr<Step>& get_current_step() const = 0;
    [[nodiscard]] virtual const shared_ptr<Converger>& get_current_converger() const = 0;
    [[nodiscard]] virtual const shared_ptr<Integrator>& get_current_integrator() const = 0;
    [[nodiscard]] virtual const shared_ptr<Solver>& get_current_solver() const = 0;

    /**
     * \brief concurrently safe insertion method
     */
    virtual void insert_loaded_dof(const uvec&) = 0;
    /**
     * \brief concurrently safe insertion method
     */
    virtual void insert_restrained_dof(const uvec&) = 0;
    /**
     * \brief concurrently safe insertion method
     */
    virtual void insert_constrained_dof(const uvec&) = 0;

    /**
     * \brief concurrently safe insertion method
     */
    virtual void insert_loaded_dof(uword) = 0;
    /**
     * \brief concurrently safe insertion method
     */
    virtual void insert_restrained_dof(uword) = 0;
    /**
     * \brief concurrently safe insertion method
     */
    virtual void insert_constrained_dof(uword) = 0;

    [[nodiscard]] virtual const suanpan::unordered_set<uword>& get_loaded_dof() const = 0;
    [[nodiscard]] virtual const suanpan::unordered_set<uword>& get_restrained_dof() const = 0;
    [[nodiscard]] virtual const suanpan::unordered_set<uword>& get_constrained_dof() const = 0;

    [[nodiscard]] virtual bool is_updated() const = 0;
    [[nodiscard]] virtual bool is_sparse() const = 0;

    virtual void set_color_model(ColorMethod) = 0;
    [[nodiscard]] virtual const std::vector<std::vector<unsigned>>& get_color_map() const = 0;
    [[nodiscard]] virtual std::pair<std::vector<unsigned>, suanpan::graph<unsigned>> get_element_connectivity(bool) = 0;

    virtual int reorder_dof() = 0;
    virtual int assign_color() = 0;

    virtual int restart() = 0;
    virtual int soft_restart() = 0;
    virtual int initialize() = 0;
    virtual int initialize_load() = 0;
    virtual int initialize_constraint() = 0;
    virtual int initialize_reference() = 0;

    [[nodiscard]] virtual int process_load(bool) = 0;
    [[nodiscard]] virtual int process_constraint(bool) = 0;
    [[nodiscard]] virtual int process_criterion() = 0;
    [[nodiscard]] virtual int process_modifier() = 0;

    virtual void record() = 0;
    virtual void enable_all() = 0;
    virtual void summary() const = 0;

    virtual void update_current_resistance() const = 0;
    virtual void update_current_damping_force() const = 0;
    virtual void update_current_inertial_force() const = 0;

    virtual void assemble_resistance() const = 0;
    virtual void assemble_damping_force() const = 0;
    virtual void assemble_inertial_force() const = 0;

    virtual void assemble_initial_mass() const = 0;
    virtual void assemble_current_mass() const = 0;
    virtual void assemble_trial_mass() const = 0;
    virtual void assemble_initial_damping() const = 0;
    virtual void assemble_current_damping() const = 0;
    virtual void assemble_trial_damping() const = 0;
    virtual void assemble_initial_stiffness() const = 0;
    virtual void assemble_current_stiffness() const = 0;
    virtual void assemble_trial_stiffness() const = 0;
    virtual void assemble_initial_geometry() const = 0;
    virtual void assemble_current_geometry() const = 0;
    virtual void assemble_trial_geometry() const = 0;

    virtual void assemble_mass_container() const = 0;
    virtual void assemble_stiffness_container() const = 0;

    virtual void erase_machine_error(vec&) const = 0;

    virtual void update_load() = 0;
    virtual void update_constraint() = 0;

    virtual void assemble_load_stiffness() = 0;
    virtual void assemble_constraint_stiffness() = 0;

    [[nodiscard]] virtual int update_current_status() const = 0;
    [[nodiscard]] virtual int update_incre_status() const = 0;
    [[nodiscard]] virtual int update_trial_status() const = 0;

    virtual void stage_status() = 0;
    virtual void commit_status() const = 0;
    virtual void clear_status() = 0;
    virtual void reset_status() const = 0;

    virtual void save(string) = 0;
};

template<typename T> bool DomainBase::erase(unsigned) { throw invalid_argument("unsupported"); }

template<> inline bool DomainBase::erase<Amplitude>(const unsigned T) { return erase_amplitude(T); }

template<> inline bool DomainBase::erase<Constraint>(const unsigned T) { return erase_constraint(T); }

template<> inline bool DomainBase::erase<Converger>(const unsigned T) { return erase_converger(T); }

template<> inline bool DomainBase::erase<Criterion>(const unsigned T) { return erase_criterion(T); }

template<> inline bool DomainBase::erase<Database>(const unsigned T) { return erase_database(T); }

template<> inline bool DomainBase::erase<Element>(const unsigned T) { return erase_element(T); }

template<> inline bool DomainBase::erase<Group>(const unsigned T) { return erase_group(T); }

template<> inline bool DomainBase::erase<Integrator>(const unsigned T) { return erase_integrator(T); }

template<> inline bool DomainBase::erase<Load>(const unsigned T) { return erase_load(T); }

template<> inline bool DomainBase::erase<Material>(const unsigned T) { return erase_material(T); }

template<> inline bool DomainBase::erase<Modifier>(const unsigned T) { return erase_modifier(T); }

template<> inline bool DomainBase::erase<Node>(const unsigned T) { return erase_node(T); }

template<> inline bool DomainBase::erase<Orientation>(const unsigned T) { return erase_orientation(T); }

template<> inline bool DomainBase::erase<Recorder>(const unsigned T) { return erase_recorder(T); }

template<> inline bool DomainBase::erase<Section>(const unsigned T) { return erase_section(T); }

template<> inline bool DomainBase::erase<Solver>(const unsigned T) { return erase_solver(T); }

template<> inline bool DomainBase::erase<Step>(const unsigned T) { return erase_step(T); }

template<typename T> const shared_ptr<T>& DomainBase::get(unsigned) { throw invalid_argument("unsupported"); }

template<typename T> const shared_ptr<T>& DomainBase::get(uword) { throw invalid_argument("unsupported"); }

template<typename T> std::vector<shared_ptr<T>> DomainBase::get(const uvec& P) {
    std::vector<shared_ptr<T>> output;
    output.reserve(P.n_elem);

    for(auto I : P) output.emplace_back(get<T>(I));

    return output;
}

template<> inline const shared_ptr<Amplitude>& DomainBase::get<Amplitude>(const uword T) { return get_amplitude(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Constraint>& DomainBase::get<Constraint>(const uword T) { return get_constraint(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Converger>& DomainBase::get<Converger>(const uword T) { return get_converger(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Criterion>& DomainBase::get<Criterion>(const uword T) { return get_criterion(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Database>& DomainBase::get<Database>(const uword T) { return get_database(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Element>& DomainBase::get<Element>(const uword T) { return get_element(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Group>& DomainBase::get<Group>(const uword T) { return get_group(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Integrator>& DomainBase::get<Integrator>(const uword T) { return get_integrator(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Load>& DomainBase::get<Load>(const uword T) { return get_load(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Material>& DomainBase::get<Material>(const uword T) { return get_material(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Modifier>& DomainBase::get<Modifier>(const uword T) { return get_modifier(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Node>& DomainBase::get<Node>(const uword T) { return get_node(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Orientation>& DomainBase::get<Orientation>(const uword T) { return get_orientation(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Recorder>& DomainBase::get<Recorder>(const uword T) { return get_recorder(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Section>& DomainBase::get<Section>(const uword T) { return get_section(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Solver>& DomainBase::get<Solver>(const uword T) { return get_solver(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Step>& DomainBase::get<Step>(const uword T) { return get_step(static_cast<unsigned>(T)); }

template<> inline const shared_ptr<Amplitude>& DomainBase::get<Amplitude>(const unsigned T) { return get_amplitude(T); }

template<> inline const shared_ptr<Constraint>& DomainBase::get<Constraint>(const unsigned T) { return get_constraint(T); }

template<> inline const shared_ptr<Converger>& DomainBase::get<Converger>(const unsigned T) { return get_converger(T); }

template<> inline const shared_ptr<Criterion>& DomainBase::get<Criterion>(const unsigned T) { return get_criterion(T); }

template<> inline const shared_ptr<Database>& DomainBase::get<Database>(const unsigned T) { return get_database(T); }

template<> inline const shared_ptr<Element>& DomainBase::get<Element>(const unsigned T) { return get_element(T); }

template<> inline const shared_ptr<Group>& DomainBase::get<Group>(const unsigned T) { return get_group(T); }

template<> inline const shared_ptr<Integrator>& DomainBase::get<Integrator>(const unsigned T) { return get_integrator(T); }

template<> inline const shared_ptr<Load>& DomainBase::get<Load>(const unsigned T) { return get_load(T); }

template<> inline const shared_ptr<Material>& DomainBase::get<Material>(const unsigned T) { return get_material(T); }

template<> inline const shared_ptr<Modifier>& DomainBase::get<Modifier>(const unsigned T) { return get_modifier(T); }

template<> inline const shared_ptr<Node>& DomainBase::get<Node>(const unsigned T) { return get_node(T); }

template<> inline const shared_ptr<Orientation>& DomainBase::get<Orientation>(const unsigned T) { return get_orientation(T); }

template<> inline const shared_ptr<Recorder>& DomainBase::get<Recorder>(const unsigned T) { return get_recorder(T); }

template<> inline const shared_ptr<Section>& DomainBase::get<Section>(const unsigned T) { return get_section(T); }

template<> inline const shared_ptr<Solver>& DomainBase::get<Solver>(const unsigned T) { return get_solver(T); }

template<> inline const shared_ptr<Step>& DomainBase::get<Step>(const unsigned T) { return get_step(T); }

template<typename T> const std::vector<shared_ptr<T>>& DomainBase::get_pool() { throw invalid_argument("unsupported"); }

template<> inline const std::vector<shared_ptr<Amplitude>>& DomainBase::get_pool<Amplitude>() { return get_amplitude_pool(); }

template<> inline const std::vector<shared_ptr<Constraint>>& DomainBase::get_pool<Constraint>() { return get_constraint_pool(); }

template<> inline const std::vector<shared_ptr<Converger>>& DomainBase::get_pool<Converger>() { return get_converger_pool(); }

template<> inline const std::vector<shared_ptr<Criterion>>& DomainBase::get_pool<Criterion>() { return get_criterion_pool(); }

template<> inline const std::vector<shared_ptr<Database>>& DomainBase::get_pool<Database>() { return get_database_pool(); }

template<> inline const std::vector<shared_ptr<Element>>& DomainBase::get_pool<Element>() { return get_element_pool(); }

template<> inline const std::vector<shared_ptr<Group>>& DomainBase::get_pool<Group>() { return get_group_pool(); }

template<> inline const std::vector<shared_ptr<Integrator>>& DomainBase::get_pool<Integrator>() { return get_integrator_pool(); }

template<> inline const std::vector<shared_ptr<Load>>& DomainBase::get_pool<Load>() { return get_load_pool(); }

template<> inline const std::vector<shared_ptr<Material>>& DomainBase::get_pool<Material>() { return get_material_pool(); }

template<> inline const std::vector<shared_ptr<Modifier>>& DomainBase::get_pool<Modifier>() { return get_modifier_pool(); }

template<> inline const std::vector<shared_ptr<Node>>& DomainBase::get_pool<Node>() { return get_node_pool(); }

template<> inline const std::vector<shared_ptr<Orientation>>& DomainBase::get_pool<Orientation>() { return get_orientation_pool(); }

template<> inline const std::vector<shared_ptr<Recorder>>& DomainBase::get_pool<Recorder>() { return get_recorder_pool(); }

template<> inline const std::vector<shared_ptr<Section>>& DomainBase::get_pool<Section>() { return get_section_pool(); }

template<> inline const std::vector<shared_ptr<Solver>>& DomainBase::get_pool<Solver>() { return get_solver_pool(); }

template<typename T> size_t DomainBase::get() { throw invalid_argument("unsupported"); }

template<> inline size_t DomainBase::get<Amplitude>() { return get_amplitude(); }

template<> inline size_t DomainBase::get<Constraint>() { return get_constraint(); }

template<> inline size_t DomainBase::get<Converger>() { return get_converger(); }

template<> inline size_t DomainBase::get<Criterion>() { return get_criterion(); }

template<> inline size_t DomainBase::get<Database>() { return get_database(); }

template<> inline size_t DomainBase::get<Element>() { return get_element(); }

template<> inline size_t DomainBase::get<Group>() { return get_group(); }

template<> inline size_t DomainBase::get<Integrator>() { return get_integrator(); }

template<> inline size_t DomainBase::get<Load>() { return get_load(); }

template<> inline size_t DomainBase::get<Material>() { return get_material(); }

template<> inline size_t DomainBase::get<Modifier>() { return get_modifier(); }

template<> inline size_t DomainBase::get<Node>() { return get_node(); }

template<> inline size_t DomainBase::get<Orientation>() { return get_orientation(); }

template<> inline size_t DomainBase::get<Recorder>() { return get_recorder(); }

template<> inline size_t DomainBase::get<Section>() { return get_section(); }

template<> inline size_t DomainBase::get<Solver>() { return get_solver(); }

template<> inline size_t DomainBase::get<Step>() { return get_step(); }

template<typename T> bool DomainBase::find(unsigned) { throw invalid_argument("unsupported"); }

template<typename T> bool DomainBase::find(uword) { throw invalid_argument("unsupported"); }

template<typename T> bool DomainBase::find(const uvec& P) {
    for(auto I : P) if(!find<T>(I)) return false;

    return true;
}

template<> inline bool DomainBase::find<Amplitude>(const uword T) { return find_amplitude(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Constraint>(const uword T) { return find_constraint(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Converger>(const uword T) { return find_converger(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Criterion>(const uword T) { return find_criterion(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Database>(const uword T) { return find_database(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Element>(const uword T) { return find_element(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Group>(const uword T) { return find_group(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Integrator>(const uword T) { return find_integrator(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Load>(const uword T) { return find_load(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Material>(const uword T) { return find_material(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Modifier>(const uword T) { return find_modifier(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Node>(const uword T) { return find_node(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Orientation>(const uword T) { return find_orientation(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Recorder>(const uword T) { return find_recorder(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Section>(const uword T) { return find_section(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Solver>(const uword T) { return find_solver(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Step>(const uword T) { return find_step(static_cast<unsigned>(T)); }

template<> inline bool DomainBase::find<Amplitude>(const unsigned T) { return find_amplitude(T); }

template<> inline bool DomainBase::find<Constraint>(const unsigned T) { return find_constraint(T); }

template<> inline bool DomainBase::find<Converger>(const unsigned T) { return find_converger(T); }

template<> inline bool DomainBase::find<Criterion>(const unsigned T) { return find_criterion(T); }

template<> inline bool DomainBase::find<Database>(const unsigned T) { return find_database(T); }

template<> inline bool DomainBase::find<Element>(const unsigned T) { return find_element(T); }

template<> inline bool DomainBase::find<Group>(const unsigned T) { return find_group(T); }

template<> inline bool DomainBase::find<Integrator>(const unsigned T) { return find_integrator(T); }

template<> inline bool DomainBase::find<Load>(const unsigned T) { return find_load(T); }

template<> inline bool DomainBase::find<Material>(const unsigned T) { return find_material(T); }

template<> inline bool DomainBase::find<Modifier>(const unsigned T) { return find_modifier(T); }

template<> inline bool DomainBase::find<Node>(const unsigned T) { return find_node(T); }

template<> inline bool DomainBase::find<Orientation>(const unsigned T) { return find_orientation(T); }

template<> inline bool DomainBase::find<Recorder>(const unsigned T) { return find_recorder(T); }

template<> inline bool DomainBase::find<Section>(const unsigned T) { return find_section(T); }

template<> inline bool DomainBase::find<Solver>(const unsigned T) { return find_solver(T); }

template<> inline bool DomainBase::find<Step>(const unsigned T) { return find_step(T); }

#endif

//! @}
