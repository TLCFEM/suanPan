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
/**
 * @class Domain
 * @brief A Domain class holds all FE model components.
 * @author tlc
 * @date 01/10/2017
 * @version 0.3.1
 * @file Domain.h
 * @addtogroup Domain
 * @{
 */

#ifndef DOMAIN_H
#define DOMAIN_H

#include <Domain/DomainBase.h>
#include <Domain/Storage.hpp>
#include <array>

using ExternalModuleQueue = std::vector<shared_ptr<ExternalModule>>;
using ThreadQueue = std::vector<std::future<void>>;

class Domain final : public DomainBase, public std::enable_shared_from_this<Domain> {
    std::atomic_bool updated = false;
    ColorMethod color_model = ColorMethod::MIS;

    unsigned current_step_tag = 0;
    std::pair<unsigned, unsigned> current_converger_tag{0, 0};  // current converger tag, current step tag
    std::pair<unsigned, unsigned> current_integrator_tag{0, 0}; // current integrator tag, current step tag
    std::pair<unsigned, unsigned> current_solver_tag{0, 0};     // current solver tag, current step tag

    ThreadQueue thread_pond;

    // dynamic libraries should be destroyed after all dependent objects are destroyed
    ExternalModuleQueue external_module_pond;

    shared_ptr<LongFactory> factory; /**< working room */

    StepQueue step_pond;

    AmplitudeStorage amplitude_pond;
    ExpressionStorage expression_pond;
    ConstraintStorage constraint_pond;
    ConvergerStorage converger_pond;
    CriterionStorage criterion_pond;
    DatabaseStorage database_pond;
    ElementStorage element_pond;
    GroupStorage group_pond;
    IntegratorStorage integrator_pond;
    LoadStorage load_pond;
    MaterialStorage material_pond;
    ModifierStorage modifier_pond;
    NodeStorage node_pond;
    OrientationStorage orientation_pond;
    RecorderStorage recorder_pond;
    SectionStorage section_pond;
    SolverStorage solver_pond;

    suanpan::unordered_set<uword> constrained_dofs; /**< data storage */
    suanpan::unordered_set<uword> loaded_dofs;      /**< data storage */
    suanpan::unordered_set<uword> restrained_dofs;  /**< data storage */

    std::vector<std::vector<unsigned>> color_map;

    std::vector<bool> attribute;

    mutable std::array<double, 5> statistics{};

public:
    explicit Domain(unsigned = 0);
    Domain(const Domain&) = delete;            // copy forbidden
    Domain(Domain&&) = delete;                 // move forbidden
    Domain& operator=(const Domain&) = delete; // assign forbidden
    Domain& operator=(Domain&&) = delete;      // assign forbidden
    ~Domain() override;

    void set_factory(const shared_ptr<LongFactory>&) override;
    const shared_ptr<LongFactory>& get_factory() const override;

    bool insert(std::future<void>&&) override;

    void wait() override;

    bool insert(const shared_ptr<ExternalModule>&) override;
    const ExternalModuleQueue& get_external_module_pool() const override;

    bool insert(const shared_ptr<Amplitude>&) override;
    bool insert(const shared_ptr<Expression>&) override;
    bool insert(const shared_ptr<Constraint>&) override;
    bool insert(const shared_ptr<Converger>&) override;
    bool insert(const shared_ptr<Criterion>&) override;
    bool insert(const shared_ptr<Database>&) override;
    bool insert(const shared_ptr<Element>&) override;
    bool insert(const shared_ptr<Group>&) override;
    bool insert(const shared_ptr<Integrator>&) override;
    bool insert(const shared_ptr<Load>&) override;
    bool insert(const shared_ptr<Material>&) override;
    bool insert(const shared_ptr<Modifier>&) override;
    bool insert(const shared_ptr<Node>&) override;
    bool insert(const shared_ptr<Orientation>&) override;
    bool insert(const shared_ptr<Recorder>&) override;
    bool insert(const shared_ptr<Section>&) override;
    bool insert(const shared_ptr<Solver>&) override;
    bool insert(const shared_ptr<Step>&) override;

    bool erase_amplitude(unsigned) override;
    bool erase_expression(unsigned) override;
    bool erase_constraint(unsigned) override;
    bool erase_converger(unsigned) override;
    bool erase_criterion(unsigned) override;
    bool erase_database(unsigned) override;
    bool erase_element(unsigned) override;
    bool erase_group(unsigned) override;
    bool erase_integrator(unsigned) override;
    bool erase_load(unsigned) override;
    bool erase_material(unsigned) override;
    bool erase_modifier(unsigned) override;
    bool erase_node(unsigned) override;
    bool erase_orientation(unsigned) override;
    bool erase_recorder(unsigned) override;
    bool erase_section(unsigned) override;
    bool erase_solver(unsigned) override;
    bool erase_step(unsigned) override;

    void disable_amplitude(unsigned) override;
    void disable_expression(unsigned) override;
    void disable_constraint(unsigned) override;
    void disable_converger(unsigned) override;
    void disable_criterion(unsigned) override;
    void disable_database(unsigned) override;
    void disable_element(unsigned) override;
    void disable_group(unsigned) override;
    void disable_integrator(unsigned) override;
    void disable_load(unsigned) override;
    void disable_material(unsigned) override;
    void disable_modifier(unsigned) override;
    void disable_node(unsigned) override;
    void disable_orientation(unsigned) override;
    void disable_recorder(unsigned) override;
    void disable_section(unsigned) override;
    void disable_solver(unsigned) override;
    void disable_step(unsigned) override;

    void enable_amplitude(unsigned) override;
    void enable_expression(unsigned) override;
    void enable_constraint(unsigned) override;
    void enable_converger(unsigned) override;
    void enable_criterion(unsigned) override;
    void enable_database(unsigned) override;
    void enable_element(unsigned) override;
    void enable_group(unsigned) override;
    void enable_integrator(unsigned) override;
    void enable_load(unsigned) override;
    void enable_material(unsigned) override;
    void enable_modifier(unsigned) override;
    void enable_node(unsigned) override;
    void enable_orientation(unsigned) override;
    void enable_recorder(unsigned) override;
    void enable_section(unsigned) override;
    void enable_solver(unsigned) override;
    void enable_step(unsigned) override;

    const shared_ptr<Amplitude>& get_amplitude(unsigned) const override;
    const shared_ptr<Expression>& get_expression(unsigned) const override;
    const shared_ptr<Constraint>& get_constraint(unsigned) const override;
    const shared_ptr<Converger>& get_converger(unsigned) const override;
    const shared_ptr<Criterion>& get_criterion(unsigned) const override;
    const shared_ptr<Database>& get_database(unsigned) const override;
    const shared_ptr<Element>& get_element(unsigned) const override;
    const shared_ptr<Group>& get_group(unsigned) const override;
    const shared_ptr<Integrator>& get_integrator(unsigned) const override;
    const shared_ptr<Load>& get_load(unsigned) const override;
    const shared_ptr<Material>& get_material(unsigned) const override;
    const shared_ptr<Modifier>& get_modifier(unsigned) const override;
    const shared_ptr<Node>& get_node(unsigned) const override;
    const shared_ptr<Orientation>& get_orientation(unsigned) const override;
    const shared_ptr<Recorder>& get_recorder(unsigned) const override;
    const shared_ptr<Section>& get_section(unsigned) const override;
    const shared_ptr<Solver>& get_solver(unsigned) const override;
    const shared_ptr<Step>& get_step(unsigned) const override;

    const AmplitudeQueue& get_amplitude_pool() const override;
    const ExpressionQueue& get_expression_pool() const override;
    const ConstraintQueue& get_constraint_pool() const override;
    const ConvergerQueue& get_converger_pool() const override;
    const CriterionQueue& get_criterion_pool() const override;
    const DatabaseQueue& get_database_pool() const override;
    const ElementQueue& get_element_pool() const override;
    const GroupQueue& get_group_pool() const override;
    const IntegratorQueue& get_integrator_pool() const override;
    const LoadQueue& get_load_pool() const override;
    const MaterialQueue& get_material_pool() const override;
    const ModifierQueue& get_modifier_pool() const override;
    const NodeQueue& get_node_pool() const override;
    const OrientationQueue& get_orientation_pool() const override;
    const RecorderQueue& get_recorder_pool() const override;
    const SectionQueue& get_section_pool() const override;
    const SolverQueue& get_solver_pool() const override;
    const StepQueue& get_step_pool() const override;

    friend shared_ptr<Amplitude>& get_amplitude(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Expression>& get_expression(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Constraint>& get_constraint(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Converger>& get_converger(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Criterion>& get_criterion(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Database>& get_database(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Element>& get_element(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Group>& get_group(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Integrator>& get_integrator(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Load>& get_load(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Material>& get_material(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Modifier>& get_modifier(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Node>& get_node(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Orientation>& get_orientation(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Recorder>& get_recorder(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Section>& get_section(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Solver>& get_solver(const shared_ptr<Domain>&, unsigned);
    friend shared_ptr<Step>& get_step(const shared_ptr<Domain>&, unsigned);

    friend shared_ptr<Amplitude>& get_amplitude(const shared_ptr<DomainBase>&, unsigned);
    friend shared_ptr<Expression>& get_expression(const shared_ptr<DomainBase>&, unsigned);
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

    size_t get_amplitude() const override;
    size_t get_expression() const override;
    size_t get_constraint() const override;
    size_t get_converger() const override;
    size_t get_criterion() const override;
    size_t get_database() const override;
    size_t get_element() const override;
    size_t get_group() const override;
    size_t get_integrator() const override;
    size_t get_load() const override;
    size_t get_material() const override;
    size_t get_modifier() const override;
    size_t get_node() const override;
    size_t get_orientation() const override;
    size_t get_recorder() const override;
    size_t get_section() const override;
    size_t get_solver() const override;
    size_t get_step() const override;

    bool find_amplitude(unsigned) const override;
    bool find_expression(unsigned) const override;
    bool find_constraint(unsigned) const override;
    bool find_converger(unsigned) const override;
    bool find_criterion(unsigned) const override;
    bool find_database(unsigned) const override;
    bool find_element(unsigned) const override;
    bool find_group(unsigned) const override;
    bool find_integrator(unsigned) const override;
    bool find_load(unsigned) const override;
    bool find_material(unsigned) const override;
    bool find_modifier(unsigned) const override;
    bool find_node(unsigned) const override;
    bool find_orientation(unsigned) const override;
    bool find_recorder(unsigned) const override;
    bool find_section(unsigned) const override;
    bool find_solver(unsigned) const override;
    bool find_step(unsigned) const override;

    void set_current_step_tag(unsigned) override;
    void set_current_converger_tag(unsigned) override;
    void set_current_integrator_tag(unsigned) override;
    void set_current_solver_tag(unsigned) override;

    unsigned get_current_step_tag() override;
    std::pair<unsigned, unsigned> get_current_converger_tag() override;
    std::pair<unsigned, unsigned> get_current_integrator_tag() override;
    std::pair<unsigned, unsigned> get_current_solver_tag() override;

    const shared_ptr<Step>& get_current_step() const override;
    const shared_ptr<Converger>& get_current_converger() const override;
    const shared_ptr<Integrator>& get_current_integrator() const override;
    const shared_ptr<Solver>& get_current_solver() const override;

    unique_ptr<Amplitude> initialized_amplitude_copy(uword) override;
    unique_ptr<Material> initialized_material_copy(uword) override;
    unique_ptr<Section> initialized_section_copy(uword) override;

    void insert_loaded_dof(const uvec&) override;
    void insert_restrained_dof(const uvec&) override;
    void insert_constrained_dof(const uvec&) override;

    void insert_loaded_dof(uword) override;
    void insert_restrained_dof(uword) override;
    void insert_constrained_dof(uword) override;

    const suanpan::unordered_set<uword>& get_loaded_dof() const override;
    const suanpan::unordered_set<uword>& get_restrained_dof() const override;
    const suanpan::unordered_set<uword>& get_constrained_dof() const override;

    bool is_updated() const override;
    bool is_sparse() const override;

    void set_attribute(ModalAttribute) override;
    [[nodiscard]] bool get_attribute(ModalAttribute) override;

    void set_color_model(ColorMethod) override;
    const std::vector<std::vector<unsigned>>& get_color_map() const override;
    std::pair<std::vector<unsigned>, suanpan::graph<unsigned>> get_element_connectivity(bool) override;

    int reorder_dof() override;
    int assign_color() override;

    // restart domain from the previous step
    int restart() override;
    // initialize domain from the previous step
    int soft_restart() override;
    // initialize the domain
    int initialize() override;
    // initialize loads for each step
    int initialize_load() override;
    // initialize constraints for each step
    int initialize_constraint() override;
    // initialize constraints for each step
    int initialize_reference() override;
    // initialize materials for each step
    int initialize_material() override;
    // initialize sections for each step
    int initialize_section() override;
    // process loads and constraints
    [[nodiscard]] int process_load(bool) override;
    [[nodiscard]] int process_constraint(bool) override;
    [[nodiscard]] int process_criterion() override;
    [[nodiscard]] int process_modifier() override;
    // record response
    void record() override;
    // enable all objects
    void enable_all() override;
    // print out domain summary
    void summary() const override;

    void update_current_resistance() const override;
    void update_current_damping_force() const override;
    void update_current_nonviscous_force() const override;
    void update_current_inertial_force() const override;

    void assemble_resistance() const override;
    void assemble_damping_force() const override;
    void assemble_nonviscous_force() const override;
    void assemble_inertial_force() const override;

    void assemble_initial_mass() const override;
    void assemble_current_mass() const override;
    void assemble_trial_mass() const override;
    void assemble_initial_damping() const override;
    void assemble_current_damping() const override;
    void assemble_trial_damping() const override;
    void assemble_initial_nonviscous() const override;
    void assemble_current_nonviscous() const override;
    void assemble_trial_nonviscous() const override;
    void assemble_initial_stiffness() const override;
    void assemble_current_stiffness() const override;
    void assemble_trial_stiffness() const override;
    void assemble_initial_geometry() const override;
    void assemble_current_geometry() const override;
    void assemble_trial_geometry() const override;

    void assemble_mass_container() const override;
    void assemble_stiffness_container() const override;

    void erase_machine_error(vec&) const override;

    void update_load() override;
    void update_constraint() override;

    void assemble_load_stiffness() override;
    void assemble_constraint_stiffness() override;

    int update_current_status() const override;
    int update_incre_status() const override;
    int update_trial_status() const override;

    void stage_status() override;
    void commit_status() const override;
    void clear_status() override;
    void reset_status() const override;

    void update(const Statistics T, const double value) const override { statistics[static_cast<size_t>(T)] += value; }

    double stats(const Statistics T) const override { return statistics[static_cast<size_t>(T)]; }

    void save(std::string) override;
};

#endif

//! @}
