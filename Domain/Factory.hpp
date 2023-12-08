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
/**
 * @class Factory
 * @brief A Factory class.
 *
 * The Factory class.
 *
 * @author tlc
 * @date 13/09/2017
 * @version 0.1.0
 * @file Factory.hpp
 * @addtogroup Storage
 * @{
 */

#ifndef FACTORY_HPP
#define FACTORY_HPP

#include <future>
#include <Toolbox/container.h>
#include <Element/MappingDOF.h>
#include <Domain/MetaMat/MetaMat>

#ifdef SUANPAN_MAGMA
#include <magmasparse.h>
#endif

enum class AnalysisType {
    NONE,
    DISP,
    EIGEN,
    BUCKLE,
    STATICS,
    DYNAMICS
};

enum class StorageScheme {
    FULL,
    BAND,
    BANDSYMM,
    SYMMPACK,
    SPARSE,
    SPARSESYMM
};

enum class SolverType {
    LAPACK,
    SPIKE,
    SUPERLU,
    MUMPS,
    CUDA,
    PARDISO,
    FGMRES,
    MAGMA,
    LIS
};

template<sp_d T> class Factory final {
    unsigned n_size = 0;               // number of degrees of freedom
    unsigned n_lobw = 0;               // low bandwidth
    unsigned n_upbw = 0;               // up bandwidth
    unsigned n_sfbw = n_lobw + n_upbw; // matrix storage offset
    unsigned n_rfld = 0;               // reference load size
    unsigned n_mpc = 0;                // multipoint constraint size
    uword n_elem = 0;

    AnalysisType analysis_type = AnalysisType::NONE;  // type of analysis
    StorageScheme storage_type = StorageScheme::FULL; // type of analysis

#ifdef SUANPAN_MAGMA
    magma_dopts magma_setting{};
#endif

    bool nlgeom = false;
    bool nonviscous = false;

    SolverType solver = SolverType::LAPACK;
    SolverSetting<T> setting{};

    T error = T(0); // error produced by certain solvers

    Col<T> ninja; // the result from A*X=B
    Col<T> sushi; // modified right-hand side B

    suanpan::set<uword> reference_dof;
    SpMat<T> reference_load;

    uvec auxiliary_encoding;      // for constraints using multiplier method
    Col<T> auxiliary_lambda;      // for constraints using multiplier method
    Col<T> auxiliary_resistance;  // for constraints using multiplier method
    Col<T> auxiliary_load;        // for constraints using multiplier method
    SpMat<T> auxiliary_stiffness; // for constraints using multiplier method

    SpCol<T> trial_constraint_resistance;
    SpCol<T> current_constraint_resistance;

    T trial_time = T(0);   // global trial (pseudo) time
    T incre_time = T(0);   // global incremental (pseudo) time
    T current_time = T(0); // global current (pseudo) time
    T pre_time = T(0);     // global previous (pseudo) time

    T strain_energy = T(0);
    T kinetic_energy = T(0);
    T viscous_energy = T(0);
    T nonviscous_energy = T(0);
    T complementary_energy = T(0);
    Col<T> momentum;

    Col<T> trial_load_factor;      // global trial load factor
    Col<T> trial_load;             // global trial load vector
    Col<T> trial_settlement;       // global trial displacement load vector
    Col<T> trial_resistance;       // global trial resistance vector
    Col<T> trial_damping_force;    // global trial damping force vector
    Col<T> trial_nonviscous_force; // global trial nonviscous damping force vector
    Col<T> trial_inertial_force;   // global trial inertial force vector
    Col<T> trial_displacement;     // global trial displacement vector
    Col<T> trial_velocity;         // global trial velocity vector
    Col<T> trial_acceleration;     // global trial acceleration vector
    Col<T> trial_temperature;      // global trial temperature vector

    Col<T> incre_load_factor;      // global incremental load vector
    Col<T> incre_load;             // global incremental load vector
    Col<T> incre_settlement;       // global incremental displacement load vector
    Col<T> incre_resistance;       // global incremental resistance vector
    Col<T> incre_damping_force;    // global incremental damping force vector
    Col<T> incre_nonviscous_force; // global incremental nonviscous damping force vector
    Col<T> incre_inertial_force;   // global incremental inertial force vector
    Col<T> incre_displacement;     // global incremental displacement vector
    Col<T> incre_velocity;         // global incremental velocity vector
    Col<T> incre_acceleration;     // global incremental acceleration vector
    Col<T> incre_temperature;      // global incremental temperature vector

    Col<T> current_load_factor;      // global current load vector
    Col<T> current_load;             // global current load vector
    Col<T> current_settlement;       // global current displacement load vector
    Col<T> current_resistance;       // global current resistance vector
    Col<T> current_damping_force;    // global current damping force vector
    Col<T> current_nonviscous_force; // global current nonviscous damping force vector
    Col<T> current_inertial_force;   // global current inertial force vector
    Col<T> current_displacement;     // global current displacement vector
    Col<T> current_velocity;         // global current velocity vector
    Col<T> current_acceleration;     // global current acceleration vector
    Col<T> current_temperature;      // global current temperature vector

    Col<T> pre_load_factor;      // global previous load vector
    Col<T> pre_load;             // global previous load vector
    Col<T> pre_settlement;       // global previous displacement load vector
    Col<T> pre_resistance;       // global previous resistance vector
    Col<T> pre_damping_force;    // global previous damping force vector
    Col<T> pre_nonviscous_force; // global previous nonviscous damping force vector
    Col<T> pre_inertial_force;   // global previous inertial force vector
    Col<T> pre_displacement;     // global previous displacement vector
    Col<T> pre_velocity;         // global previous velocity vector
    Col<T> pre_acceleration;     // global previous acceleration vector
    Col<T> pre_temperature;      // global previous temperature vector

    shared_ptr<MetaMat<T>> global_mass = nullptr;       // global mass matrix
    shared_ptr<MetaMat<T>> global_damping = nullptr;    // global damping matrix
    shared_ptr<MetaMat<T>> global_nonviscous = nullptr; // global nonviscous damping matrix
    shared_ptr<MetaMat<T>> global_stiffness = nullptr;  // global stiffness matrix
    shared_ptr<MetaMat<T>> global_geometry = nullptr;   // global geometry matrix

    std::vector<std::mutex> global_mutex{20};

    Col<T> eigenvalue; // eigenvalues

    Mat<T> eigenvector; // eigenvectors

    unique_ptr<MetaMat<T>> get_basic_container();
    unique_ptr<MetaMat<T>> get_matrix_container();

    void assemble_matrix_helper(shared_ptr<MetaMat<T>>&, const Mat<T>&, const uvec&, const std::vector<MappingDOF>&);

public:
    const bool initialized = false;

    explicit Factory(unsigned = 0, AnalysisType = AnalysisType::NONE, StorageScheme = StorageScheme::FULL);

    void set_size(unsigned);
    [[nodiscard]] unsigned get_size() const;

    void set_entry(uword);
    [[nodiscard]] uword get_entry() const;

    void set_nlgeom(bool);
    [[nodiscard]] bool is_nlgeom() const;

    void set_nonviscous(bool);
    [[nodiscard]] bool is_nonviscous() const;

    void set_solver_type(SolverType);
    [[nodiscard]] SolverType get_solver_type() const;

    void set_solver_setting(const SolverSetting<double>&);
    [[nodiscard]] const SolverSetting<double>& get_solver_setting() const;

#ifdef SUANPAN_MAGMA
    void set_solver_setting(const magma_dopts& magma_opt) { magma_setting = magma_opt; }

    [[nodiscard]] const magma_dopts& get_magma_setting() const { return magma_setting; }
#endif

    void set_analysis_type(AnalysisType);
    [[nodiscard]] AnalysisType get_analysis_type() const;

    void set_storage_scheme(StorageScheme);
    [[nodiscard]] StorageScheme get_storage_scheme() const;

    [[nodiscard]] bool is_sparse() const;

    void set_bandwidth(unsigned, unsigned);
    void get_bandwidth(unsigned&, unsigned&) const;

    void update_reference_size();
    void set_reference_size(unsigned);
    [[nodiscard]] unsigned get_reference_size() const;

    void update_reference_dof(const uvec&);
    void set_reference_dof(const suanpan::set<uword>&);
    [[nodiscard]] const suanpan::set<uword>& get_reference_dof() const;

    void set_error(T);
    T get_error() const;

    /*************************INITIALIZER*************************/

    int initialize();

    void initialize_load_factor();
    void initialize_load();
    void initialize_settlement();
    void initialize_resistance();
    void initialize_damping_force();
    void initialize_nonviscous_force();
    void initialize_inertial_force();
    void initialize_displacement();
    void initialize_velocity();
    void initialize_acceleration();
    void initialize_temperature();
    void initialize_auxiliary_resistance();

    void initialize_mass();
    void initialize_damping();
    void initialize_nonviscous();
    void initialize_stiffness();
    void initialize_geometry();
    void initialize_eigen();

    /*************************SETTER*************************/

    void set_ninja(const Col<T>&);
    void set_sushi(const Col<T>&);

    void update_sushi_by(const Col<T>&);

    void set_mpc(unsigned);

    void set_reference_load(const SpMat<T>&);

    void set_trial_time(T);
    void set_trial_load_factor(const Col<T>&);
    void set_trial_load(const Col<T>&);
    void set_trial_settlement(const Col<T>&);
    void set_trial_resistance(const Col<T>&);
    void set_trial_damping_force(const Col<T>&);
    void set_trial_nonviscous_force(const Col<T>&);
    void set_trial_inertial_force(const Col<T>&);
    void set_trial_displacement(const Col<T>&);
    void set_trial_velocity(const Col<T>&);
    void set_trial_acceleration(const Col<T>&);
    void set_trial_temperature(const Col<T>&);

    void set_incre_time(T);
    void set_incre_load_factor(const Col<T>&);
    void set_incre_load(const Col<T>&);
    void set_incre_settlement(const Col<T>&);
    void set_incre_resistance(const Col<T>&);
    void set_incre_damping_force(const Col<T>&);
    void set_incre_nonviscous_force(const Col<T>&);
    void set_incre_inertial_force(const Col<T>&);
    void set_incre_displacement(const Col<T>&);
    void set_incre_velocity(const Col<T>&);
    void set_incre_acceleration(const Col<T>&);
    void set_incre_temperature(const Col<T>&);

    void set_current_time(T);
    void set_current_load_factor(const Col<T>&);
    void set_current_load(const Col<T>&);
    void set_current_settlement(const Col<T>&);
    void set_current_resistance(const Col<T>&);
    void set_current_damping_force(const Col<T>&);
    void set_current_nonviscous_force(const Col<T>&);
    void set_current_inertial_force(const Col<T>&);
    void set_current_displacement(const Col<T>&);
    void set_current_velocity(const Col<T>&);
    void set_current_acceleration(const Col<T>&);
    void set_current_temperature(const Col<T>&);

    void set_pre_time(T);
    void set_pre_load_factor(const Col<T>&);
    void set_pre_load(const Col<T>&);
    void set_pre_settlement(const Col<T>&);
    void set_pre_resistance(const Col<T>&);
    void set_pre_damping_force(const Col<T>&);
    void set_pre_nonviscous_force(const Col<T>&);
    void set_pre_inertial_force(const Col<T>&);
    void set_pre_displacement(const Col<T>&);
    void set_pre_velocity(const Col<T>&);
    void set_pre_acceleration(const Col<T>&);
    void set_pre_temperature(const Col<T>&);

    void set_mass(const shared_ptr<MetaMat<T>>&);
    void set_damping(const shared_ptr<MetaMat<T>>&);
    void set_nonviscous(const shared_ptr<MetaMat<T>>&);
    void set_stiffness(const shared_ptr<MetaMat<T>>&);
    void set_geometry(const shared_ptr<MetaMat<T>>&);

    void set_eigenvalue(const Col<T>&);
    void set_eigenvector(const Mat<T>&);

    /*************************GETTER*************************/

    const Col<T>& get_ninja() const;
    const Col<T>& get_sushi() const;

    [[nodiscard]] unsigned get_mpc() const;

    const SpMat<T>& get_reference_load() const;

    [[nodiscard]] const uvec& get_auxiliary_encoding() const;
    const Col<T>& get_auxiliary_lambda() const;
    const Col<T>& get_auxiliary_resistance() const;
    const Col<T>& get_auxiliary_load() const;
    const SpMat<T>& get_auxiliary_stiffness() const;

    const SpCol<T>& get_trial_constraint_resistance() const;
    const SpCol<T>& get_current_constraint_resistance() const;

    T get_strain_energy();
    T get_kinetic_energy();
    T get_viscous_energy();
    T get_nonviscous_energy();
    T get_complementary_energy();
    const Col<T>& get_momentum();

    T get_trial_time() const;
    const Col<T>& get_trial_load_factor() const;
    const Col<T>& get_trial_load() const;
    const Col<T>& get_trial_settlement() const;
    const Col<T>& get_trial_resistance() const;
    const Col<T>& get_trial_damping_force() const;
    const Col<T>& get_trial_nonviscous_force() const;
    const Col<T>& get_trial_inertial_force() const;
    const Col<T>& get_trial_displacement() const;
    const Col<T>& get_trial_velocity() const;
    const Col<T>& get_trial_acceleration() const;
    const Col<T>& get_trial_temperature() const;

    T get_incre_time() const;
    const Col<T>& get_incre_load_factor() const;
    const Col<T>& get_incre_load() const;
    const Col<T>& get_incre_settlement() const;
    const Col<T>& get_incre_resistance() const;
    const Col<T>& get_incre_damping_force() const;
    const Col<T>& get_incre_nonviscous_force() const;
    const Col<T>& get_incre_inertial_force() const;
    const Col<T>& get_incre_displacement() const;
    const Col<T>& get_incre_velocity() const;
    const Col<T>& get_incre_acceleration() const;
    const Col<T>& get_incre_temperature() const;

    T get_current_time() const;
    const Col<T>& get_current_load_factor() const;
    const Col<T>& get_current_load() const;
    const Col<T>& get_current_settlement() const;
    const Col<T>& get_current_resistance() const;
    const Col<T>& get_current_damping_force() const;
    const Col<T>& get_current_nonviscous_force() const;
    const Col<T>& get_current_inertial_force() const;
    const Col<T>& get_current_displacement() const;
    const Col<T>& get_current_velocity() const;
    const Col<T>& get_current_acceleration() const;
    const Col<T>& get_current_temperature() const;

    T get_pre_time() const;
    const Col<T>& get_pre_load_factor() const;
    const Col<T>& get_pre_load() const;
    const Col<T>& get_pre_settlement() const;
    const Col<T>& get_pre_resistance() const;
    const Col<T>& get_pre_damping_force() const;
    const Col<T>& get_pre_nonviscous_force() const;
    const Col<T>& get_pre_inertial_force() const;
    const Col<T>& get_pre_displacement() const;
    const Col<T>& get_pre_velocity() const;
    const Col<T>& get_pre_acceleration() const;
    const Col<T>& get_pre_temperature() const;

    const shared_ptr<MetaMat<T>>& get_mass() const;
    const shared_ptr<MetaMat<T>>& get_damping() const;
    const shared_ptr<MetaMat<T>>& get_nonviscous() const;
    const shared_ptr<MetaMat<T>>& get_stiffness() const;
    const shared_ptr<MetaMat<T>>& get_geometry() const;

    std::mutex& get_auxiliary_encoding_mutex();
    std::mutex& get_auxiliary_resistance_mutex();
    std::mutex& get_auxiliary_load_mutex();
    std::mutex& get_auxiliary_stiffness_mutex();

    std::mutex& get_trial_constraint_resistance_mutex();

    std::mutex& get_trial_load_mutex();
    std::mutex& get_trial_settlement_mutex();
    std::mutex& get_reference_load_mutex();

    std::mutex& get_mass_mutex();
    std::mutex& get_damping_mutex();
    std::mutex& get_nonviscous_mutex();
    std::mutex& get_stiffness_mutex();
    std::mutex& get_geometry_mutex();

    const Col<T>& get_eigenvalue() const;
    const Mat<T>& get_eigenvector() const;

    /*************************UPDATER*************************/

    void update_trial_time(T);
    void update_trial_load_factor(const Col<T>&);
    void update_trial_load(const Col<T>&);
    void update_trial_settlement(const Col<T>&);
    void update_trial_resistance(const Col<T>&);
    void update_trial_damping_force(const Col<T>&);
    void update_trial_nonviscous_force(const Col<T>&);
    void update_trial_inertial_force(const Col<T>&);
    void update_trial_displacement(const Col<T>&);
    void update_trial_velocity(const Col<T>&);
    void update_trial_acceleration(const Col<T>&);
    void update_trial_temperature(const Col<T>&);

    void update_incre_time(T);
    void update_incre_load_factor(const Col<T>&);
    void update_incre_load(const Col<T>&);
    void update_incre_settlement(const Col<T>&);
    void update_incre_resistance(const Col<T>&);
    void update_incre_damping_force(const Col<T>&);
    void update_incre_nonviscous_force(const Col<T>&);
    void update_incre_inertial_force(const Col<T>&);
    void update_incre_displacement(const Col<T>&);
    void update_incre_velocity(const Col<T>&);
    void update_incre_acceleration(const Col<T>&);
    void update_incre_temperature(const Col<T>&);

    void update_current_time(T);
    void update_current_load_factor(const Col<T>&);
    void update_current_load(const Col<T>&);
    void update_current_settlement(const Col<T>&);
    void update_current_resistance(const Col<T>&);
    void update_current_damping_force(const Col<T>&);
    void update_current_nonviscous_force(const Col<T>&);
    void update_current_inertial_force(const Col<T>&);
    void update_current_displacement(const Col<T>&);
    void update_current_velocity(const Col<T>&);
    void update_current_acceleration(const Col<T>&);
    void update_current_temperature(const Col<T>&);

    void update_trial_time_by(T);
    void update_trial_load_factor_by(const Col<T>&);
    void update_trial_load_by(const Col<T>&);
    void update_trial_settlement_by(const Col<T>&);
    void update_trial_resistance_by(const Col<T>&);
    void update_trial_damping_force_by(const Col<T>&);
    void update_trial_nonviscous_force_by(const Col<T>&);
    void update_trial_inertial_force_by(const Col<T>&);
    void update_trial_displacement_by(const Col<T>&);
    void update_trial_velocity_by(const Col<T>&);
    void update_trial_acceleration_by(const Col<T>&);
    void update_trial_temperature_by(const Col<T>&);

    void update_incre_time_by(T);
    void update_incre_load_factor_by(const Col<T>&);
    void update_incre_load_by(const Col<T>&);
    void update_incre_settlement_by(const Col<T>&);
    void update_incre_resistance_by(const Col<T>&);
    void update_incre_damping_force_by(const Col<T>&);
    void update_incre_nonviscous_force_by(const Col<T>&);
    void update_incre_inertial_force_by(const Col<T>&);
    void update_incre_displacement_by(const Col<T>&);
    void update_incre_velocity_by(const Col<T>&);
    void update_incre_acceleration_by(const Col<T>&);
    void update_incre_temperature_by(const Col<T>&);

    void update_current_time_by(T);
    void update_current_load_factor_by(const Col<T>&);
    void update_current_load_by(const Col<T>&);
    void update_current_settlement_by(const Col<T>&);
    void update_current_resistance_by(const Col<T>&);
    void update_current_damping_force_by(const Col<T>&);
    void update_current_nonviscous_force_by(const Col<T>&);
    void update_current_inertial_force_by(const Col<T>&);
    void update_current_displacement_by(const Col<T>&);
    void update_current_velocity_by(const Col<T>&);
    void update_current_acceleration_by(const Col<T>&);
    void update_current_temperature_by(const Col<T>&);

    /*************************FRIEND*************************/

    Col<T>& modify_ninja();
    Col<T>& modify_sushi();

    suanpan::set<uword>& modify_reference_dof();
    SpMat<T>& modify_reference_load();

    uvec& modify_auxiliary_encoding();
    Col<T>& modify_auxiliary_lambda();
    Col<T>& modify_auxiliary_resistance();
    Col<T>& modify_auxiliary_load();
    SpMat<T>& modify_auxiliary_stiffness();

    SpCol<T>& modify_trial_constraint_resistance();
    SpCol<T>& modify_current_constraint_resistance();

    T& modify_trial_time();
    Col<T>& modify_trial_load_factor();
    Col<T>& modify_trial_load();
    Col<T>& modify_trial_settlement();
    Col<T>& modify_trial_resistance();
    Col<T>& modify_trial_damping_force();
    Col<T>& modify_trial_nonviscous_force();
    Col<T>& modify_trial_inertial_force();
    Col<T>& modify_trial_displacement();
    Col<T>& modify_trial_velocity();
    Col<T>& modify_trial_acceleration();
    Col<T>& modify_trial_temperature();

    T& modify_incre_time();
    Col<T>& modify_incre_load_factor();
    Col<T>& modify_incre_load();
    Col<T>& modify_incre_settlement();
    Col<T>& modify_incre_resistance();
    Col<T>& modify_incre_damping_force();
    Col<T>& modify_incre_nonviscous_force();
    Col<T>& modify_incre_inertial_force();
    Col<T>& modify_incre_displacement();
    Col<T>& modify_incre_velocity();
    Col<T>& modify_incre_acceleration();
    Col<T>& modify_incre_temperature();

    T& modify_current_time();
    Col<T>& modify_current_load_factor();
    Col<T>& modify_current_load();
    Col<T>& modify_current_settlement();
    Col<T>& modify_current_resistance();
    Col<T>& modify_current_damping_force();
    Col<T>& modify_current_nonviscous_force();
    Col<T>& modify_current_inertial_force();
    Col<T>& modify_current_displacement();
    Col<T>& modify_current_velocity();
    Col<T>& modify_current_acceleration();
    Col<T>& modify_current_temperature();

    T& modify_pre_time();
    Col<T>& modify_pre_load_factor();
    Col<T>& modify_pre_load();
    Col<T>& modify_pre_settlement();
    Col<T>& modify_pre_resistance();
    Col<T>& modify_pre_damping_force();
    Col<T>& modify_pre_nonviscous_force();
    Col<T>& modify_pre_inertial_force();
    Col<T>& modify_pre_displacement();
    Col<T>& modify_pre_velocity();
    Col<T>& modify_pre_acceleration();
    Col<T>& modify_pre_temperature();

    shared_ptr<MetaMat<T>>& modify_mass();
    shared_ptr<MetaMat<T>>& modify_damping();
    shared_ptr<MetaMat<T>>& modify_nonviscous();
    shared_ptr<MetaMat<T>>& modify_stiffness();
    shared_ptr<MetaMat<T>>& modify_geometry();

    Col<T>& modify_eigenvalue();
    Mat<T>& modify_eigenvector();

    /*************************STATUS*************************/

    void commit_energy();
    void clear_energy();

    void commit_status();
    void commit_time();
    void commit_load_factor();
    void commit_load();
    void commit_settlement();
    void commit_resistance();
    void commit_damping_force();
    void commit_nonviscous_force();
    void commit_inertial_force();
    void commit_displacement();
    void commit_velocity();
    void commit_acceleration();
    void commit_temperature();
    void commit_auxiliary_resistance();

    void commit_pre_status();
    void commit_pre_time();
    void commit_pre_load_factor();
    void commit_pre_load();
    void commit_pre_settlement();
    void commit_pre_resistance();
    void commit_pre_damping_force();
    void commit_pre_nonviscous_force();
    void commit_pre_inertial_force();
    void commit_pre_displacement();
    void commit_pre_velocity();
    void commit_pre_acceleration();
    void commit_pre_temperature();

    void clear_status();
    void clear_time();
    void clear_load_factor();
    void clear_load();
    void clear_settlement();
    void clear_resistance();
    void clear_damping_force();
    void clear_nonviscous_force();
    void clear_inertial_force();
    void clear_displacement();
    void clear_velocity();
    void clear_acceleration();
    void clear_temperature();
    void clear_auxiliary_resistance();

    void reset_status();
    void reset_time();
    void reset_load_factor();
    void reset_load();
    void reset_settlement();
    void reset_resistance();
    void reset_damping_force();
    void reset_nonviscous_force();
    void reset_inertial_force();
    void reset_displacement();
    void reset_velocity();
    void reset_acceleration();
    void reset_temperature();
    void reset_auxiliary_resistance();

    void clear_eigen();
    void clear_mass();
    void clear_damping();
    void clear_nonviscous();
    void clear_stiffness();
    void clear_geometry();
    void clear_auxiliary();

    void reset();

    /*************************ASSEMBLER*************************/

    void assemble_resistance(const Mat<T>&, const uvec&);
    void assemble_damping_force(const Mat<T>&, const uvec&);
    void assemble_nonviscous_force(const Mat<T>&, const uvec&);
    void assemble_inertial_force(const Mat<T>&, const uvec&);

    void assemble_mass(const Mat<T>&, const uvec&, const std::vector<MappingDOF>&);
    void assemble_damping(const Mat<T>&, const uvec&, const std::vector<MappingDOF>&);
    void assemble_nonviscous(const Mat<T>&, const uvec&, const std::vector<MappingDOF>&);
    void assemble_stiffness(const Mat<T>&, const uvec&, const std::vector<MappingDOF>&);
    void assemble_geometry(const Mat<T>&, const uvec&, const std::vector<MappingDOF>&);

    void assemble_stiffness(const SpMat<T>&, const uvec&);

    /*************************UTILITY*************************/

    void print() const;
};

template<sp_d T> Factory<T>::Factory(const unsigned D, const AnalysisType AT, const StorageScheme SS)
    : n_size(D)
    , analysis_type(AT)
    , storage_type(SS) {}

template<sp_d T> void Factory<T>::set_size(const unsigned D) {
    if(D == n_size) return;
    n_size = D;
    access::rw(initialized) = false;
}

template<sp_d T> unsigned Factory<T>::get_size() const { return n_size; }

template<sp_d T> void Factory<T>::set_entry(const uword N) {
    n_elem = N;
    if(n_elem > std::numeric_limits<int>::max()) throw invalid_argument("too many elements");
}

template<sp_d T> uword Factory<T>::get_entry() const { return n_elem; }

template<sp_d T> void Factory<T>::set_nlgeom(const bool B) {
    if(B == nlgeom) return;
    nlgeom = B;
    access::rw(initialized) = false;
}

template<sp_d T> bool Factory<T>::is_nlgeom() const { return nlgeom; }

template<sp_d T> void Factory<T>::set_nonviscous(const bool B) {
    if(B == nonviscous) return;
    nonviscous = B;
    access::rw(initialized) = false;
}

template<sp_d T> bool Factory<T>::is_nonviscous() const { return nonviscous; }

template<sp_d T> void Factory<T>::set_solver_type(const SolverType E) { solver = E; }

template<sp_d T> SolverType Factory<T>::get_solver_type() const { return solver; }

template<sp_d T> void Factory<T>::set_solver_setting(const SolverSetting<double>& SS) { setting = SS; }

template<sp_d T> const SolverSetting<double>& Factory<T>::get_solver_setting() const { return setting; }

template<sp_d T> void Factory<T>::set_analysis_type(const AnalysisType AT) {
    if(AT == analysis_type) return;
    analysis_type = AT;
    access::rw(initialized) = false;
}

template<sp_d T> AnalysisType Factory<T>::get_analysis_type() const { return analysis_type; }

template<sp_d T> void Factory<T>::set_storage_scheme(const StorageScheme SS) {
    if(SS == storage_type) return;
    storage_type = SS;
    access::rw(initialized) = false;
}

template<sp_d T> StorageScheme Factory<T>::get_storage_scheme() const { return storage_type; }

template<sp_d T> bool Factory<T>::is_sparse() const { return StorageScheme::SPARSE == storage_type || StorageScheme::SPARSESYMM == storage_type; }

template<sp_d T> void Factory<T>::set_bandwidth(const unsigned L, const unsigned U) {
    if(L == n_lobw && U == n_upbw) return;
    n_lobw = L;
    n_upbw = U;
    n_sfbw = L + U;
    access::rw(initialized) = false;
}

template<sp_d T> void Factory<T>::get_bandwidth(unsigned& L, unsigned& U) const {
    L = n_lobw;
    U = n_upbw;
}

template<sp_d T> void Factory<T>::update_reference_size() { n_rfld = static_cast<unsigned>(reference_dof.size()); }

template<sp_d T> void Factory<T>::set_reference_size(const unsigned S) {
    if(S == n_rfld) return;
    n_rfld = S;
}

template<sp_d T> unsigned Factory<T>::get_reference_size() const { return n_rfld; }

template<sp_d T> void Factory<T>::update_reference_dof(const uvec& S) { reference_dof.insert(S.cbegin(), S.cend()); }

template<sp_d T> void Factory<T>::set_reference_dof(const suanpan::set<uword>& D) { reference_dof = D; }

template<sp_d T> const suanpan::set<uword>& Factory<T>::get_reference_dof() const { return reference_dof; }

template<sp_d T> void Factory<T>::set_error(const T E) { error = E; }

template<sp_d T> T Factory<T>::get_error() const { return error; }

template<sp_d T> int Factory<T>::initialize() {
    reference_dof.clear(); // clear reference dof vector in every step

    if(initialized || n_size == 0) return 0;

    ninja.zeros(n_size);
    sushi.zeros(n_size);

    reset();

    switch(analysis_type) {
    case AnalysisType::DISP:
        initialize_displacement();
        break;
    case AnalysisType::EIGEN:
        initialize_mass();
        initialize_stiffness();
        initialize_eigen();
        break;
    case AnalysisType::BUCKLE:
    case AnalysisType::STATICS:
        initialize_load();
        initialize_resistance();
        initialize_displacement();
        initialize_stiffness();
        initialize_geometry();
        break;
    case AnalysisType::DYNAMICS:
        initialize_load();
        initialize_resistance();
        initialize_damping_force();
        initialize_nonviscous_force();
        initialize_inertial_force();
        initialize_displacement();
        initialize_velocity();
        initialize_acceleration();
        initialize_mass();
        initialize_damping();
        initialize_nonviscous();
        initialize_stiffness();
        initialize_geometry();
        break;
    case AnalysisType::NONE:
        break;
    }

    initialize_auxiliary_resistance();

    access::rw(initialized) = true;

    return 0;
}

template<sp_d T> void Factory<T>::initialize_load_factor() {
    if(n_rfld == 0) return;

    trial_load_factor.zeros(n_rfld);
    incre_load_factor.zeros(n_rfld);
    current_load_factor.zeros(n_rfld);

    reference_load.zeros(n_size, n_rfld);
}

template<sp_d T> void Factory<T>::initialize_load() {
    trial_load.zeros(n_size);
    incre_load.zeros(n_size);
    current_load.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_settlement() {
    trial_settlement.zeros(n_size);
    incre_settlement.zeros(n_size);
    current_settlement.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_resistance() {
    trial_resistance.zeros(n_size);
    incre_resistance.zeros(n_size);
    current_resistance.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_damping_force() {
    trial_damping_force.zeros(n_size);
    incre_damping_force.zeros(n_size);
    current_damping_force.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_nonviscous_force() {
    trial_nonviscous_force.zeros(n_size);
    incre_nonviscous_force.zeros(n_size);
    current_nonviscous_force.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_inertial_force() {
    trial_inertial_force.zeros(n_size);
    incre_inertial_force.zeros(n_size);
    current_inertial_force.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_displacement() {
    trial_displacement.zeros(n_size);
    incre_displacement.zeros(n_size);
    current_displacement.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_velocity() {
    trial_velocity.zeros(n_size);
    incre_velocity.zeros(n_size);
    current_velocity.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_acceleration() {
    trial_acceleration.zeros(n_size);
    incre_acceleration.zeros(n_size);
    current_acceleration.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_temperature() {
    trial_temperature.zeros(n_size);
    incre_temperature.zeros(n_size);
    current_temperature.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_auxiliary_resistance() {
    trial_constraint_resistance.zeros(n_size);
    current_constraint_resistance.zeros(n_size);
}

template<sp_d T> void Factory<T>::initialize_mass() { global_mass = get_matrix_container(); }

template<sp_d T> void Factory<T>::initialize_damping() { global_damping = get_matrix_container(); }

template<sp_d T> void Factory<T>::initialize_nonviscous() {
    if(!nonviscous) return;

    global_nonviscous = get_matrix_container();
}

template<sp_d T> void Factory<T>::initialize_stiffness() { global_stiffness = get_matrix_container(); }

template<sp_d T> void Factory<T>::initialize_geometry() {
    if(!nlgeom) return;

    global_geometry = get_matrix_container();
}

template<sp_d T> void Factory<T>::initialize_eigen() {
    eigenvalue.zeros(n_size);
    eigenvector.zeros(n_size, n_size);
}

template<sp_d T> void Factory<T>::set_ninja(const Col<T>& N) { ninja = N; }

template<sp_d T> void Factory<T>::set_sushi(const Col<T>& S) { sushi = S; }

template<sp_d T> void Factory<T>::update_sushi_by(const Col<T>& S) { sushi += S; }

template<sp_d T> void Factory<T>::set_mpc(const unsigned S) {
    n_mpc = S;
    auxiliary_encoding.zeros(n_mpc);
    auxiliary_resistance.zeros(n_mpc);
    auxiliary_load.zeros(n_mpc);
    auxiliary_stiffness.zeros(n_size, n_mpc);
}

template<sp_d T> void Factory<T>::set_reference_load(const SpMat<T>& L) { reference_load = L; }

template<sp_d T> void Factory<T>::set_mass(const shared_ptr<MetaMat<T>>& M) { global_mass = M; }

template<sp_d T> void Factory<T>::set_damping(const shared_ptr<MetaMat<T>>& C) { global_damping = C; }

template<sp_d T> void Factory<T>::set_nonviscous(const shared_ptr<MetaMat<T>>& C) { global_nonviscous = C; }

template<sp_d T> void Factory<T>::set_stiffness(const shared_ptr<MetaMat<T>>& K) { global_stiffness = K; }

template<sp_d T> void Factory<T>::set_geometry(const shared_ptr<MetaMat<T>>& G) { global_geometry = G; }

template<sp_d T> void Factory<T>::set_eigenvalue(const Col<T>& L) { eigenvalue = L; }

template<sp_d T> void Factory<T>::set_eigenvector(const Mat<T>& V) { eigenvector = V; }

template<sp_d T> const Col<T>& Factory<T>::get_ninja() const { return ninja; }

template<sp_d T> const Col<T>& Factory<T>::get_sushi() const { return sushi; }

template<sp_d T> unsigned Factory<T>::get_mpc() const { return n_mpc; }

template<sp_d T> const SpMat<T>& Factory<T>::get_reference_load() const { return reference_load; }

template<sp_d T> const uvec& Factory<T>::get_auxiliary_encoding() const { return auxiliary_encoding; }

template<sp_d T> const Col<T>& Factory<T>::get_auxiliary_lambda() const { return auxiliary_lambda; }

template<sp_d T> const Col<T>& Factory<T>::get_auxiliary_resistance() const { return auxiliary_resistance; }

template<sp_d T> const Col<T>& Factory<T>::get_auxiliary_load() const { return auxiliary_load; }

template<sp_d T> const SpMat<T>& Factory<T>::get_auxiliary_stiffness() const { return auxiliary_stiffness; }

template<sp_d T> const SpCol<T>& Factory<T>::get_trial_constraint_resistance() const { return trial_constraint_resistance; }

template<sp_d T> const SpCol<T>& Factory<T>::get_current_constraint_resistance() const { return current_constraint_resistance; }

template<sp_d T> T Factory<T>::get_strain_energy() { return strain_energy; }

template<sp_d T> T Factory<T>::get_kinetic_energy() { return kinetic_energy; }

template<sp_d T> T Factory<T>::get_viscous_energy() { return viscous_energy; }

template<sp_d T> T Factory<T>::get_nonviscous_energy() { return nonviscous_energy; }

template<sp_d T> T Factory<T>::get_complementary_energy() { return complementary_energy; }

template<sp_d T> const Col<T>& Factory<T>::get_momentum() { return momentum; }

template<sp_d T> const shared_ptr<MetaMat<T>>& Factory<T>::get_mass() const { return global_mass; }

template<sp_d T> const shared_ptr<MetaMat<T>>& Factory<T>::get_damping() const { return global_damping; }

template<sp_d T> const shared_ptr<MetaMat<T>>& Factory<T>::get_nonviscous() const { return global_nonviscous; }

template<sp_d T> const shared_ptr<MetaMat<T>>& Factory<T>::get_stiffness() const { return global_stiffness; }

template<sp_d T> const shared_ptr<MetaMat<T>>& Factory<T>::get_geometry() const { return global_geometry; }

template<sp_d T> std::mutex& Factory<T>::get_auxiliary_encoding_mutex() { return global_mutex.at(0); }

template<sp_d T> std::mutex& Factory<T>::get_auxiliary_resistance_mutex() { return global_mutex.at(1); }

template<sp_d T> std::mutex& Factory<T>::get_auxiliary_load_mutex() { return global_mutex.at(2); }

template<sp_d T> std::mutex& Factory<T>::get_auxiliary_stiffness_mutex() { return global_mutex.at(3); }

template<sp_d T> std::mutex& Factory<T>::get_trial_constraint_resistance_mutex() { return global_mutex.at(4); }

template<sp_d T> std::mutex& Factory<T>::get_trial_load_mutex() { return global_mutex.at(5); }

template<sp_d T> std::mutex& Factory<T>::get_trial_settlement_mutex() { return global_mutex.at(6); }

template<sp_d T> std::mutex& Factory<T>::get_reference_load_mutex() { return global_mutex.at(7); }

template<sp_d T> std::mutex& Factory<T>::get_mass_mutex() { return global_mutex.at(8); }

template<sp_d T> std::mutex& Factory<T>::get_damping_mutex() { return global_mutex.at(9); }

template<sp_d T> std::mutex& Factory<T>::get_nonviscous_mutex() { return global_mutex.at(10); }

template<sp_d T> std::mutex& Factory<T>::get_stiffness_mutex() { return global_mutex.at(11); }

template<sp_d T> std::mutex& Factory<T>::get_geometry_mutex() { return global_mutex.at(12); }

template<sp_d T> const Col<T>& Factory<T>::get_eigenvalue() const { return eigenvalue; }

template<sp_d T> const Mat<T>& Factory<T>::get_eigenvector() const { return eigenvector; }

template<sp_d T> void Factory<T>::commit_energy() {
    auto se = std::async([&] { if(!trial_resistance.empty() && !incre_displacement.empty()) strain_energy += .5 * dot(trial_resistance + current_resistance, incre_displacement); });
    auto ke = std::async([&] { if(!trial_inertial_force.empty() && !trial_velocity.empty()) kinetic_energy = .5 * dot(global_mass * trial_velocity, trial_velocity); });
    auto ve = std::async([&] { if(!trial_damping_force.empty() && !incre_displacement.empty()) viscous_energy += .5 * dot(trial_damping_force + current_damping_force, incre_displacement); });
    auto ne = std::async([&] { if(!trial_nonviscous_force.empty() && !incre_displacement.empty()) nonviscous_energy += .5 * dot(trial_nonviscous_force + current_nonviscous_force, incre_displacement); });
    auto ce = std::async([&] { if(!trial_displacement.empty() && !incre_resistance.empty()) complementary_energy += .5 * dot(trial_displacement + current_displacement, incre_resistance); });
    auto mm = std::async([&] { if(!trial_inertial_force.empty() && !trial_velocity.empty()) momentum = global_mass * trial_velocity; });

    se.get();
    ke.get();
    ve.get();
    ne.get();
    ce.get();
    mm.get();
}

template<sp_d T> void Factory<T>::clear_energy() {
    strain_energy = T(0);
    kinetic_energy = T(0);
    viscous_energy = T(0);
    nonviscous_energy = T(0);
    complementary_energy = T(0);
    momentum.zeros();
}

template<sp_d T> void Factory<T>::commit_status() {
    ninja.zeros();

    commit_energy();

    commit_time();
    commit_load_factor();
    commit_load();
    commit_settlement();
    commit_resistance();
    commit_damping_force();
    commit_nonviscous_force();
    commit_inertial_force();
    commit_displacement();
    commit_velocity();
    commit_acceleration();
    commit_temperature();
    commit_auxiliary_resistance();
}

template<sp_d T> void Factory<T>::commit_time() {
    current_time = trial_time;
    incre_time = T(0);
}

template<sp_d T> void Factory<T>::commit_load_factor() {
    if(trial_load_factor.is_empty()) return;
    current_load_factor = trial_load_factor;
    incre_load_factor.zeros();
}

template<sp_d T> void Factory<T>::commit_load() {
    if(trial_load.is_empty()) return;
    current_load = trial_load;
    incre_load.zeros();
}

template<sp_d T> void Factory<T>::commit_settlement() {
    if(trial_settlement.is_empty()) return;
    current_settlement = trial_settlement;
    incre_settlement.zeros();
}

template<sp_d T> void Factory<T>::commit_resistance() {
    if(trial_resistance.is_empty()) return;
    current_resistance = trial_resistance;
    incre_resistance.zeros();
}

template<sp_d T> void Factory<T>::commit_damping_force() {
    if(trial_damping_force.is_empty()) return;
    current_damping_force = trial_damping_force;
    incre_damping_force.zeros();
}

template<sp_d T> void Factory<T>::commit_nonviscous_force() {
    if(trial_nonviscous_force.is_empty()) return;
    current_nonviscous_force = trial_nonviscous_force;
    incre_nonviscous_force.zeros();
}

template<sp_d T> void Factory<T>::commit_inertial_force() {
    if(trial_inertial_force.is_empty()) return;
    current_inertial_force = trial_inertial_force;
    incre_inertial_force.zeros();
}

template<sp_d T> void Factory<T>::commit_displacement() {
    if(trial_displacement.is_empty()) return;
    current_displacement = trial_displacement;
    incre_displacement.zeros();
}

template<sp_d T> void Factory<T>::commit_velocity() {
    if(trial_velocity.is_empty()) return;
    current_velocity = trial_velocity;
    incre_velocity.zeros();
}

template<sp_d T> void Factory<T>::commit_acceleration() {
    if(trial_acceleration.is_empty()) return;
    current_acceleration = trial_acceleration;
    incre_acceleration.zeros();
}

template<sp_d T> void Factory<T>::commit_temperature() {
    if(trial_temperature.is_empty()) return;
    current_temperature = trial_temperature;
    incre_temperature.zeros();
}

template<sp_d T> void Factory<T>::commit_auxiliary_resistance() {
    if(trial_constraint_resistance.is_empty()) return;
    current_constraint_resistance = trial_constraint_resistance;
}

template<sp_d T> void Factory<T>::commit_pre_status() {
    commit_pre_time();
    commit_pre_load_factor();
    commit_pre_load();
    commit_pre_settlement();
    commit_pre_resistance();
    commit_pre_damping_force();
    commit_pre_nonviscous_force();
    commit_pre_inertial_force();
    commit_pre_displacement();
    commit_pre_velocity();
    commit_pre_acceleration();
    commit_pre_temperature();
}

template<sp_d T> void Factory<T>::commit_pre_time() { pre_time = current_time; }

template<sp_d T> void Factory<T>::commit_pre_load_factor() { if(!current_load_factor.is_empty()) pre_load_factor = current_load_factor; }

template<sp_d T> void Factory<T>::commit_pre_load() { if(!current_load.is_empty()) pre_load = current_load; }

template<sp_d T> void Factory<T>::commit_pre_settlement() { if(!current_settlement.is_empty()) pre_settlement = current_settlement; }

template<sp_d T> void Factory<T>::commit_pre_resistance() { if(!current_resistance.is_empty()) pre_resistance = current_resistance; }

template<sp_d T> void Factory<T>::commit_pre_damping_force() { if(!current_damping_force.is_empty()) pre_damping_force = current_damping_force; }

template<sp_d T> void Factory<T>::commit_pre_nonviscous_force() { if(!current_nonviscous_force.is_empty()) pre_nonviscous_force = current_nonviscous_force; }

template<sp_d T> void Factory<T>::commit_pre_inertial_force() { if(!current_inertial_force.is_empty()) pre_inertial_force = current_inertial_force; }

template<sp_d T> void Factory<T>::commit_pre_displacement() { if(!current_displacement.is_empty()) pre_displacement = current_displacement; }

template<sp_d T> void Factory<T>::commit_pre_velocity() { if(!current_velocity.is_empty()) pre_velocity = current_velocity; }

template<sp_d T> void Factory<T>::commit_pre_acceleration() { if(!current_acceleration.is_empty()) pre_acceleration = current_acceleration; }

template<sp_d T> void Factory<T>::commit_pre_temperature() { if(!current_temperature.is_empty()) pre_temperature = current_temperature; }

template<sp_d T> void Factory<T>::clear_status() {
    access::rw(initialized) = false;

    ninja.zeros();

    clear_energy();

    clear_time();
    clear_load_factor();
    clear_load();
    clear_settlement();
    clear_resistance();
    clear_damping_force();
    clear_nonviscous_force();
    clear_inertial_force();
    clear_displacement();
    clear_velocity();
    clear_acceleration();
    clear_temperature();
    clear_auxiliary_resistance();
}

template<sp_d T> void Factory<T>::clear_time() { trial_time = incre_time = current_time = 0.; }

template<sp_d T> void Factory<T>::clear_load_factor() {
    if(!pre_load_factor.is_empty()) pre_load_factor.zeros();
    if(!trial_load_factor.is_empty()) trial_load_factor.zeros();
    if(!incre_load_factor.is_empty()) incre_load_factor.zeros();
    if(!current_load_factor.is_empty()) current_load_factor.zeros();
}

template<sp_d T> void Factory<T>::clear_load() {
    if(!pre_load.is_empty()) pre_load.zeros();
    if(!trial_load.is_empty()) trial_load.zeros();
    if(!incre_load.is_empty()) incre_load.zeros();
    if(!current_load.is_empty()) current_load.zeros();
}

template<sp_d T> void Factory<T>::clear_settlement() {
    if(!pre_settlement.is_empty()) pre_settlement.zeros();
    if(!trial_settlement.is_empty()) trial_settlement.zeros();
    if(!incre_settlement.is_empty()) incre_settlement.zeros();
    if(!current_settlement.is_empty()) current_settlement.zeros();
}

template<sp_d T> void Factory<T>::clear_resistance() {
    if(!pre_resistance.is_empty()) pre_resistance.zeros();
    if(!trial_resistance.is_empty()) trial_resistance.zeros();
    if(!incre_resistance.is_empty()) incre_resistance.zeros();
    if(!current_resistance.is_empty()) current_resistance.zeros();
}

template<sp_d T> void Factory<T>::clear_damping_force() {
    if(!pre_damping_force.is_empty()) pre_damping_force.zeros();
    if(!trial_damping_force.is_empty()) trial_damping_force.zeros();
    if(!incre_damping_force.is_empty()) incre_damping_force.zeros();
    if(!current_damping_force.is_empty()) current_damping_force.zeros();
}

template<sp_d T> void Factory<T>::clear_nonviscous_force() {
    if(!pre_nonviscous_force.is_empty()) pre_nonviscous_force.zeros();
    if(!trial_nonviscous_force.is_empty()) trial_nonviscous_force.zeros();
    if(!incre_nonviscous_force.is_empty()) incre_nonviscous_force.zeros();
    if(!current_nonviscous_force.is_empty()) current_nonviscous_force.zeros();
}

template<sp_d T> void Factory<T>::clear_inertial_force() {
    if(!pre_inertial_force.is_empty()) pre_inertial_force.zeros();
    if(!trial_inertial_force.is_empty()) trial_inertial_force.zeros();
    if(!incre_inertial_force.is_empty()) incre_inertial_force.zeros();
    if(!current_inertial_force.is_empty()) current_inertial_force.zeros();
}

template<sp_d T> void Factory<T>::clear_displacement() {
    if(!pre_displacement.is_empty()) pre_displacement.zeros();
    if(!trial_displacement.is_empty()) trial_displacement.zeros();
    if(!incre_displacement.is_empty()) incre_displacement.zeros();
    if(!current_displacement.is_empty()) current_displacement.zeros();
}

template<sp_d T> void Factory<T>::clear_velocity() {
    if(!pre_velocity.is_empty()) pre_velocity.zeros();
    if(!trial_velocity.is_empty()) trial_velocity.zeros();
    if(!incre_velocity.is_empty()) incre_velocity.zeros();
    if(!current_velocity.is_empty()) current_velocity.zeros();
}

template<sp_d T> void Factory<T>::clear_acceleration() {
    if(!pre_acceleration.is_empty()) pre_acceleration.zeros();
    if(!trial_acceleration.is_empty()) trial_acceleration.zeros();
    if(!incre_acceleration.is_empty()) incre_acceleration.zeros();
    if(!current_acceleration.is_empty()) current_acceleration.zeros();
}

template<sp_d T> void Factory<T>::clear_temperature() {
    if(!pre_temperature.is_empty()) pre_temperature.zeros();
    if(!trial_temperature.is_empty()) trial_temperature.zeros();
    if(!incre_temperature.is_empty()) incre_temperature.zeros();
    if(!current_temperature.is_empty()) current_temperature.zeros();
}

template<sp_d T> void Factory<T>::clear_auxiliary_resistance() {
    if(!trial_constraint_resistance.is_empty()) trial_constraint_resistance.zeros();
    if(!current_constraint_resistance.is_empty()) current_constraint_resistance.zeros();
}

template<sp_d T> void Factory<T>::reset_status() {
    ninja.zeros();

    reset_time();
    reset_load_factor();
    reset_load();
    reset_settlement();
    reset_resistance();
    reset_damping_force();
    reset_nonviscous_force();
    reset_inertial_force();
    reset_displacement();
    reset_velocity();
    reset_acceleration();
    reset_temperature();
    reset_auxiliary_resistance();
}

template<sp_d T> void Factory<T>::reset_time() {
    trial_time = current_time;
    incre_time = T(0);
}

template<sp_d T> void Factory<T>::reset_load_factor() {
    if(trial_load_factor.is_empty()) return;
    trial_load_factor = current_load_factor;
    incre_load_factor.zeros();
}

template<sp_d T> void Factory<T>::reset_load() {
    if(trial_load.is_empty()) return;
    trial_load = current_load;
    incre_load.zeros();
}

template<sp_d T> void Factory<T>::reset_settlement() {
    if(trial_settlement.is_empty()) return;
    trial_settlement = current_settlement;
    incre_settlement.zeros();
}

template<sp_d T> void Factory<T>::reset_resistance() {
    if(trial_resistance.is_empty()) return;
    trial_resistance = current_resistance;
    incre_resistance.zeros();
}

template<sp_d T> void Factory<T>::reset_damping_force() {
    if(trial_damping_force.is_empty()) return;
    trial_damping_force = current_damping_force;
    incre_damping_force.zeros();
}

template<sp_d T> void Factory<T>::reset_nonviscous_force() {
    if(trial_nonviscous_force.is_empty()) return;
    trial_nonviscous_force = current_nonviscous_force;
    incre_nonviscous_force.zeros();
}

template<sp_d T> void Factory<T>::reset_inertial_force() {
    if(trial_inertial_force.is_empty()) return;
    trial_inertial_force = current_inertial_force;
    incre_inertial_force.zeros();
}

template<sp_d T> void Factory<T>::reset_displacement() {
    if(trial_displacement.is_empty()) return;
    trial_displacement = current_displacement;
    incre_displacement.zeros();
}

template<sp_d T> void Factory<T>::reset_velocity() {
    if(trial_velocity.is_empty()) return;
    trial_velocity = current_velocity;
    incre_velocity.zeros();
}

template<sp_d T> void Factory<T>::reset_acceleration() {
    if(trial_acceleration.is_empty()) return;
    trial_acceleration = current_acceleration;
    incre_acceleration.zeros();
}

template<sp_d T> void Factory<T>::reset_temperature() {
    if(trial_temperature.is_empty()) return;
    trial_temperature = current_temperature;
    incre_temperature.zeros();
}

template<sp_d T> void Factory<T>::reset_auxiliary_resistance() {
    if(trial_constraint_resistance.is_empty()) return;
    trial_constraint_resistance = current_constraint_resistance;
}

template<sp_d T> void Factory<T>::clear_eigen() {
    if(!eigenvalue.is_empty()) eigenvalue.zeros();
    if(!eigenvector.is_empty()) eigenvector.zeros();
}

template<sp_d T> void Factory<T>::clear_mass() { if(global_mass != nullptr) global_mass->zeros(); }

template<sp_d T> void Factory<T>::clear_damping() { if(global_damping != nullptr) global_damping->zeros(); }

template<sp_d T> void Factory<T>::clear_nonviscous() { if(global_nonviscous != nullptr) global_nonviscous->zeros(); }

template<sp_d T> void Factory<T>::clear_stiffness() { if(global_stiffness != nullptr) global_stiffness->zeros(); }

template<sp_d T> void Factory<T>::clear_geometry() { if(global_geometry != nullptr) global_geometry->zeros(); }

template<sp_d T> void Factory<T>::clear_auxiliary() {
    n_mpc = 0;
    auxiliary_load.reset();
    auxiliary_stiffness.set_size(n_size, 0);
    auxiliary_resistance.reset();
    auxiliary_encoding.reset();
}

template<sp_d T> void Factory<T>::reset() {
    global_mass = nullptr;
    global_damping = nullptr;
    global_nonviscous = nullptr;
    global_stiffness = nullptr;
    global_geometry = nullptr;
}

template<sp_d T> void Factory<T>::assemble_resistance(const Mat<T>& ER, const uvec& EI) {
    if(ER.is_empty()) return;
    for(auto I = 0llu; I < EI.n_elem; ++I) trial_resistance(EI(I)) += ER(I);
}

template<sp_d T> void Factory<T>::assemble_damping_force(const Mat<T>& ER, const uvec& EI) {
    if(ER.is_empty()) return;
    for(auto I = 0llu; I < EI.n_elem; ++I) trial_damping_force(EI(I)) += ER(I);
}

template<sp_d T> void Factory<T>::assemble_nonviscous_force(const Mat<T>& ER, const uvec& EI) {
    if(ER.is_empty()) return;
    for(auto I = 0llu; I < EI.n_elem; ++I) trial_nonviscous_force(EI(I)) += ER(I);
}

template<sp_d T> void Factory<T>::assemble_inertial_force(const Mat<T>& ER, const uvec& EI) {
    if(ER.is_empty()) return;
    for(auto I = 0llu; I < EI.n_elem; ++I) trial_inertial_force(EI(I)) += ER(I);
}

template<sp_d T> void Factory<T>::assemble_matrix_helper(shared_ptr<MetaMat<T>>& GM, const Mat<T>& EM, const uvec& EI, const std::vector<MappingDOF>& MAP) {
    if(EM.is_empty()) return;

    if(StorageScheme::BANDSYMM == storage_type || StorageScheme::SYMMPACK == storage_type) for(const auto [g_row, g_col, l_row, l_col] : MAP) GM->unsafe_at(g_row, g_col) += EM(l_row, l_col);
    else for(auto I = 0llu; I < EI.n_elem; ++I) for(auto J = 0llu; J < EI.n_elem; ++J) GM->unsafe_at(EI(J), EI(I)) += EM(J, I);
}

template<sp_d T> void Factory<T>::assemble_mass(const Mat<T>& EM, const uvec& EI, const std::vector<MappingDOF>& MAP) { this->assemble_matrix_helper(global_mass, EM, EI, MAP); }

template<sp_d T> void Factory<T>::assemble_damping(const Mat<T>& EC, const uvec& EI, const std::vector<MappingDOF>& MAP) { this->assemble_matrix_helper(global_damping, EC, EI, MAP); }

template<sp_d T> void Factory<T>::assemble_nonviscous(const Mat<T>& EC, const uvec& EI, const std::vector<MappingDOF>& MAP) { this->assemble_matrix_helper(global_nonviscous, EC, EI, MAP); }

template<sp_d T> void Factory<T>::assemble_stiffness(const Mat<T>& EK, const uvec& EI, const std::vector<MappingDOF>& MAP) { this->assemble_matrix_helper(global_stiffness, EK, EI, MAP); }

template<sp_d T> void Factory<T>::assemble_geometry(const Mat<T>& EG, const uvec& EI, const std::vector<MappingDOF>& MAP) { this->assemble_matrix_helper(global_geometry, EG, EI, MAP); }

template<sp_d T> void Factory<T>::assemble_stiffness(const SpMat<T>& EK, const uvec& EI) {
    if(EK.is_empty()) return;
    for(auto I = EK.begin(); I != EK.end(); ++I) global_stiffness->at(EI(I.row()), EI(I.col())) += *I;
}

template<sp_d T> void Factory<T>::print() const {
    suanpan_info("A Factory object with size of {}.\n", n_size);
}

template<sp_d T> unique_ptr<MetaMat<T>> Factory<T>::get_basic_container() {
    switch(storage_type) {
    case StorageScheme::FULL:
#ifdef SUANPAN_CUDA
        if(SolverType::CUDA == solver) return std::make_unique<FullMatCUDA<T>>(n_size, n_size);
#endif
        return std::make_unique<FullMat<T>>(n_size, n_size);
    case StorageScheme::BAND:
        if(SolverType::SPIKE == solver) return std::make_unique<BandMatSpike<T>>(n_size, n_lobw, n_upbw);
        return std::make_unique<BandMat<T>>(n_size, n_lobw, n_upbw);
    case StorageScheme::BANDSYMM:
        return std::make_unique<BandSymmMat<T>>(n_size, n_lobw);
    case StorageScheme::SYMMPACK:
        return std::make_unique<SymmPackMat<T>>(n_size);
    case StorageScheme::SPARSE:
        if(SolverType::MUMPS == solver) return std::make_unique<SparseMatMUMPS<T>>(n_size, n_size, n_elem);
        if(SolverType::LIS == solver) return std::make_unique<SparseMatLis<T>>(n_size, n_size, n_elem);
        if(SolverType::SUPERLU == solver) return std::make_unique<SparseMatSuperLU<T>>(n_size, n_size, n_elem);
#ifdef SUANPAN_MKL
        if(SolverType::PARDISO == solver) return std::make_unique<SparseMatPARDISO<T>>(n_size, n_size, n_elem);
        if(SolverType::FGMRES == solver) return std::make_unique<SparseMatFGMRES<T>>(n_size, n_size, n_elem);
#endif
#ifdef SUANPAN_CUDA
        if(SolverType::CUDA == solver) return std::make_unique<SparseMatCUDA<T>>(n_size, n_size, n_elem);
#ifdef SUANPAN_MAGMA
        if(SolverType::MAGMA == solver) return std::make_unique<SparseMatMAGMA<T>>(n_size, n_size, magma_setting);
#endif
#endif
        return std::make_unique<SparseMatSuperLU<T>>(n_size, n_size, n_elem);
    case StorageScheme::SPARSESYMM:
#ifdef SUANPAN_MKL
        if(SolverType::FGMRES == solver) return std::make_unique<SparseSymmMatFGMRES<T>>(n_size, n_size, n_elem);
#endif
        return std::make_unique<SparseSymmMatMUMPS<T>>(n_size, n_size, n_elem);
    default:
        throw invalid_argument("need a proper storage scheme");
    }
}

template<sp_d T> unique_ptr<MetaMat<T>> Factory<T>::get_matrix_container() {
    auto global_mat = get_basic_container();

    global_mat->set_solver_setting(setting);

    return global_mat;
}

template<sp_d T> shared_ptr<MetaMat<T>>& Factory<T>::modify_mass() { return global_mass; }

template<sp_d T> shared_ptr<MetaMat<T>>& Factory<T>::modify_damping() { return global_damping; }

template<sp_d T> shared_ptr<MetaMat<T>>& Factory<T>::modify_nonviscous() { return global_nonviscous; }

template<sp_d T> shared_ptr<MetaMat<T>>& Factory<T>::modify_stiffness() { return global_stiffness; }

template<sp_d T> shared_ptr<MetaMat<T>>& Factory<T>::modify_geometry() { return global_geometry; }

template<sp_d T> Col<T>& Factory<T>::modify_ninja() { return ninja; }

template<sp_d T> Col<T>& Factory<T>::modify_sushi() { return sushi; }

template<sp_d T> suanpan::set<uword>& Factory<T>::modify_reference_dof() { return reference_dof; }

template<sp_d T> SpMat<T>& Factory<T>::modify_reference_load() { return reference_load; }

template<sp_d T> uvec& Factory<T>::modify_auxiliary_encoding() { return auxiliary_encoding; }

template<sp_d T> Col<T>& Factory<T>::modify_auxiliary_lambda() { return auxiliary_lambda; }

template<sp_d T> Col<T>& Factory<T>::modify_auxiliary_resistance() { return auxiliary_resistance; }

template<sp_d T> Col<T>& Factory<T>::modify_auxiliary_load() { return auxiliary_load; }

template<sp_d T> SpMat<T>& Factory<T>::modify_auxiliary_stiffness() { return auxiliary_stiffness; }

template<sp_d T> SpCol<T>& Factory<T>::modify_trial_constraint_resistance() { return trial_constraint_resistance; }

template<sp_d T> SpCol<T>& Factory<T>::modify_current_constraint_resistance() { return current_constraint_resistance; }

template<sp_d T> Col<T>& Factory<T>::modify_eigenvalue() { return eigenvalue; }

template<sp_d T> Mat<T>& Factory<T>::modify_eigenvector() { return eigenvector; }

template<sp_d T> void Factory<T>::set_trial_time(const T M) { trial_time = M; }

template<sp_d T> void Factory<T>::set_trial_load_factor(const Col<T>& L) { trial_load_factor = L; }

template<sp_d T> void Factory<T>::set_trial_load(const Col<T>& L) { trial_load = L; }

template<sp_d T> void Factory<T>::set_trial_settlement(const Col<T>& S) { trial_settlement = S; }

template<sp_d T> void Factory<T>::set_trial_resistance(const Col<T>& R) { trial_resistance = R; }

template<sp_d T> void Factory<T>::set_trial_damping_force(const Col<T>& R) { trial_damping_force = R; }

template<sp_d T> void Factory<T>::set_trial_nonviscous_force(const Col<T>& R) { trial_nonviscous_force = R; }

template<sp_d T> void Factory<T>::set_trial_inertial_force(const Col<T>& R) { trial_inertial_force = R; }

template<sp_d T> void Factory<T>::set_trial_displacement(const Col<T>& D) { trial_displacement = D; }

template<sp_d T> void Factory<T>::set_trial_velocity(const Col<T>& V) { trial_velocity = V; }

template<sp_d T> void Factory<T>::set_trial_acceleration(const Col<T>& A) { trial_acceleration = A; }

template<sp_d T> void Factory<T>::set_trial_temperature(const Col<T>& M) { trial_temperature = M; }

template<sp_d T> void Factory<T>::set_incre_time(const T M) { incre_time = M; }

template<sp_d T> void Factory<T>::set_incre_load_factor(const Col<T>& L) { incre_load_factor = L; }

template<sp_d T> void Factory<T>::set_incre_load(const Col<T>& L) { incre_load = L; }

template<sp_d T> void Factory<T>::set_incre_settlement(const Col<T>& S) { incre_settlement = S; }

template<sp_d T> void Factory<T>::set_incre_resistance(const Col<T>& R) { incre_resistance = R; }

template<sp_d T> void Factory<T>::set_incre_damping_force(const Col<T>& R) { incre_damping_force = R; }

template<sp_d T> void Factory<T>::set_incre_nonviscous_force(const Col<T>& R) { incre_nonviscous_force = R; }

template<sp_d T> void Factory<T>::set_incre_inertial_force(const Col<T>& R) { incre_inertial_force = R; }

template<sp_d T> void Factory<T>::set_incre_displacement(const Col<T>& D) { incre_displacement = D; }

template<sp_d T> void Factory<T>::set_incre_velocity(const Col<T>& V) { incre_velocity = V; }

template<sp_d T> void Factory<T>::set_incre_acceleration(const Col<T>& A) { incre_acceleration = A; }

template<sp_d T> void Factory<T>::set_incre_temperature(const Col<T>& M) { incre_temperature = M; }

template<sp_d T> void Factory<T>::set_current_time(const T M) { current_time = M; }

template<sp_d T> void Factory<T>::set_current_load_factor(const Col<T>& L) { current_load_factor = L; }

template<sp_d T> void Factory<T>::set_current_load(const Col<T>& L) { current_load = L; }

template<sp_d T> void Factory<T>::set_current_settlement(const Col<T>& S) { current_settlement = S; }

template<sp_d T> void Factory<T>::set_current_resistance(const Col<T>& R) { current_resistance = R; }

template<sp_d T> void Factory<T>::set_current_damping_force(const Col<T>& R) { current_damping_force = R; }

template<sp_d T> void Factory<T>::set_current_nonviscous_force(const Col<T>& R) { current_nonviscous_force = R; }

template<sp_d T> void Factory<T>::set_current_inertial_force(const Col<T>& R) { current_inertial_force = R; }

template<sp_d T> void Factory<T>::set_current_displacement(const Col<T>& D) { current_displacement = D; }

template<sp_d T> void Factory<T>::set_current_velocity(const Col<T>& V) { current_velocity = V; }

template<sp_d T> void Factory<T>::set_current_acceleration(const Col<T>& A) { current_acceleration = A; }

template<sp_d T> void Factory<T>::set_current_temperature(const Col<T>& M) { current_temperature = M; }

template<sp_d T> void Factory<T>::set_pre_time(const T M) { pre_time = M; }

template<sp_d T> void Factory<T>::set_pre_load_factor(const Col<T>& L) { pre_load_factor = L; }

template<sp_d T> void Factory<T>::set_pre_load(const Col<T>& L) { pre_load = L; }

template<sp_d T> void Factory<T>::set_pre_settlement(const Col<T>& S) { pre_settlement = S; }

template<sp_d T> void Factory<T>::set_pre_resistance(const Col<T>& R) { pre_resistance = R; }

template<sp_d T> void Factory<T>::set_pre_damping_force(const Col<T>& R) { pre_damping_force = R; }

template<sp_d T> void Factory<T>::set_pre_nonviscous_force(const Col<T>& R) { pre_nonviscous_force = R; }

template<sp_d T> void Factory<T>::set_pre_inertial_force(const Col<T>& R) { pre_inertial_force = R; }

template<sp_d T> void Factory<T>::set_pre_displacement(const Col<T>& D) { pre_displacement = D; }

template<sp_d T> void Factory<T>::set_pre_velocity(const Col<T>& V) { pre_velocity = V; }

template<sp_d T> void Factory<T>::set_pre_acceleration(const Col<T>& A) { pre_acceleration = A; }

template<sp_d T> void Factory<T>::set_pre_temperature(const Col<T>& M) { pre_temperature = M; }

template<sp_d T> T Factory<T>::get_trial_time() const { return trial_time; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_load_factor() const { return trial_load_factor; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_load() const { return trial_load; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_settlement() const { return trial_settlement; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_resistance() const { return trial_resistance; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_damping_force() const { return trial_damping_force; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_nonviscous_force() const { return trial_nonviscous_force; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_inertial_force() const { return trial_inertial_force; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_displacement() const { return trial_displacement; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_velocity() const { return trial_velocity; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_acceleration() const { return trial_acceleration; }

template<sp_d T> const Col<T>& Factory<T>::get_trial_temperature() const { return trial_temperature; }

template<sp_d T> T Factory<T>::get_incre_time() const { return incre_time; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_load_factor() const { return incre_load_factor; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_load() const { return incre_load; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_settlement() const { return incre_settlement; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_resistance() const { return incre_resistance; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_damping_force() const { return incre_damping_force; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_nonviscous_force() const { return incre_nonviscous_force; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_inertial_force() const { return incre_inertial_force; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_displacement() const { return incre_displacement; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_velocity() const { return incre_velocity; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_acceleration() const { return incre_acceleration; }

template<sp_d T> const Col<T>& Factory<T>::get_incre_temperature() const { return incre_temperature; }

template<sp_d T> T Factory<T>::get_current_time() const { return current_time; }

template<sp_d T> const Col<T>& Factory<T>::get_current_load_factor() const { return current_load_factor; }

template<sp_d T> const Col<T>& Factory<T>::get_current_load() const { return current_load; }

template<sp_d T> const Col<T>& Factory<T>::get_current_settlement() const { return current_settlement; }

template<sp_d T> const Col<T>& Factory<T>::get_current_resistance() const { return current_resistance; }

template<sp_d T> const Col<T>& Factory<T>::get_current_damping_force() const { return current_damping_force; }

template<sp_d T> const Col<T>& Factory<T>::get_current_nonviscous_force() const { return current_nonviscous_force; }

template<sp_d T> const Col<T>& Factory<T>::get_current_inertial_force() const { return current_inertial_force; }

template<sp_d T> const Col<T>& Factory<T>::get_current_displacement() const { return current_displacement; }

template<sp_d T> const Col<T>& Factory<T>::get_current_velocity() const { return current_velocity; }

template<sp_d T> const Col<T>& Factory<T>::get_current_acceleration() const { return current_acceleration; }

template<sp_d T> const Col<T>& Factory<T>::get_current_temperature() const { return current_temperature; }

template<sp_d T> T Factory<T>::get_pre_time() const { return pre_time; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_load_factor() const { return pre_load_factor; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_load() const { return pre_load; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_settlement() const { return pre_settlement; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_resistance() const { return pre_resistance; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_damping_force() const { return pre_damping_force; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_nonviscous_force() const { return pre_nonviscous_force; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_inertial_force() const { return pre_inertial_force; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_displacement() const { return pre_displacement; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_velocity() const { return pre_velocity; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_acceleration() const { return pre_acceleration; }

template<sp_d T> const Col<T>& Factory<T>::get_pre_temperature() const { return pre_temperature; }

template<sp_d T> T& Factory<T>::modify_trial_time() { return trial_time; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_load_factor() { return trial_load_factor; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_load() { return trial_load; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_settlement() { return trial_settlement; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_resistance() { return trial_resistance; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_damping_force() { return trial_damping_force; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_nonviscous_force() { return trial_nonviscous_force; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_inertial_force() { return trial_inertial_force; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_displacement() { return trial_displacement; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_velocity() { return trial_velocity; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_acceleration() { return trial_acceleration; }

template<sp_d T> Col<T>& Factory<T>::modify_trial_temperature() { return trial_temperature; }

template<sp_d T> T& Factory<T>::modify_incre_time() { return incre_time; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_load_factor() { return incre_load_factor; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_load() { return incre_load; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_settlement() { return incre_settlement; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_resistance() { return incre_resistance; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_damping_force() { return incre_damping_force; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_nonviscous_force() { return incre_nonviscous_force; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_inertial_force() { return incre_inertial_force; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_displacement() { return incre_displacement; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_velocity() { return incre_velocity; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_acceleration() { return incre_acceleration; }

template<sp_d T> Col<T>& Factory<T>::modify_incre_temperature() { return incre_temperature; }

template<sp_d T> T& Factory<T>::modify_current_time() { return current_time; }

template<sp_d T> Col<T>& Factory<T>::modify_current_load_factor() { return current_load_factor; }

template<sp_d T> Col<T>& Factory<T>::modify_current_load() { return current_load; }

template<sp_d T> Col<T>& Factory<T>::modify_current_settlement() { return current_settlement; }

template<sp_d T> Col<T>& Factory<T>::modify_current_resistance() { return current_resistance; }

template<sp_d T> Col<T>& Factory<T>::modify_current_damping_force() { return current_damping_force; }

template<sp_d T> Col<T>& Factory<T>::modify_current_nonviscous_force() { return current_nonviscous_force; }

template<sp_d T> Col<T>& Factory<T>::modify_current_inertial_force() { return current_inertial_force; }

template<sp_d T> Col<T>& Factory<T>::modify_current_displacement() { return current_displacement; }

template<sp_d T> Col<T>& Factory<T>::modify_current_velocity() { return current_velocity; }

template<sp_d T> Col<T>& Factory<T>::modify_current_acceleration() { return current_acceleration; }

template<sp_d T> Col<T>& Factory<T>::modify_current_temperature() { return current_temperature; }

template<sp_d T> T& Factory<T>::modify_pre_time() { return pre_time; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_load_factor() { return pre_load_factor; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_load() { return pre_load; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_settlement() { return pre_settlement; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_resistance() { return pre_resistance; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_damping_force() { return pre_damping_force; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_nonviscous_force() { return pre_nonviscous_force; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_inertial_force() { return pre_inertial_force; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_displacement() { return pre_displacement; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_velocity() { return pre_velocity; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_acceleration() { return pre_acceleration; }

template<sp_d T> Col<T>& Factory<T>::modify_pre_temperature() { return pre_temperature; }

template<sp_d T> void Factory<T>::update_trial_time(const T M) {
    trial_time = M;
    incre_time = trial_time - current_time;
}

template<sp_d T> void Factory<T>::update_trial_load_factor(const Col<T>& L) {
    trial_load_factor = L;
    incre_load_factor = trial_load_factor - current_load_factor;
}

template<sp_d T> void Factory<T>::update_trial_load(const Col<T>& L) {
    trial_load = L;
    incre_load = trial_load - current_load;
}

template<sp_d T> void Factory<T>::update_trial_settlement(const Col<T>& S) {
    trial_settlement = S;
    incre_settlement = trial_settlement - current_settlement;
}

template<sp_d T> void Factory<T>::update_trial_resistance(const Col<T>& R) {
    trial_resistance = R;
    incre_resistance = trial_resistance - current_resistance;
}

template<sp_d T> void Factory<T>::update_trial_damping_force(const Col<T>& R) {
    trial_damping_force = R;
    incre_damping_force = trial_damping_force - current_damping_force;
}

template<sp_d T> void Factory<T>::update_trial_nonviscous_force(const Col<T>& R) {
    trial_nonviscous_force = R;
    incre_nonviscous_force = trial_nonviscous_force - current_nonviscous_force;
}

template<sp_d T> void Factory<T>::update_trial_inertial_force(const Col<T>& R) {
    trial_inertial_force = R;
    incre_inertial_force = trial_inertial_force - current_inertial_force;
}

template<sp_d T> void Factory<T>::update_trial_displacement(const Col<T>& D) {
    trial_displacement = D;
    incre_displacement = trial_displacement - current_displacement;
}

template<sp_d T> void Factory<T>::update_trial_velocity(const Col<T>& V) {
    trial_velocity = V;
    incre_velocity = trial_velocity - current_velocity;
}

template<sp_d T> void Factory<T>::update_trial_acceleration(const Col<T>& A) {
    trial_acceleration = A;
    incre_acceleration = trial_acceleration - current_acceleration;
}

template<sp_d T> void Factory<T>::update_trial_temperature(const Col<T>& M) {
    trial_temperature = M;
    incre_temperature = trial_temperature - current_temperature;
}

template<sp_d T> void Factory<T>::update_incre_time(const T M) {
    incre_time = M;
    trial_time = current_time + incre_time;
}

template<sp_d T> void Factory<T>::update_incre_load_factor(const Col<T>& L) {
    incre_load_factor = L;
    trial_load_factor = current_load_factor + incre_load_factor;
}

template<sp_d T> void Factory<T>::update_incre_load(const Col<T>& L) {
    incre_load = L;
    trial_load = current_load + incre_load;
}

template<sp_d T> void Factory<T>::update_incre_settlement(const Col<T>& S) {
    incre_settlement = S;
    trial_settlement = current_settlement + incre_settlement;
}

template<sp_d T> void Factory<T>::update_incre_resistance(const Col<T>& R) {
    incre_resistance = R;
    trial_resistance = current_resistance + incre_resistance;
}

template<sp_d T> void Factory<T>::update_incre_damping_force(const Col<T>& R) {
    incre_damping_force = R;
    trial_damping_force = current_damping_force + incre_damping_force;
}

template<sp_d T> void Factory<T>::update_incre_nonviscous_force(const Col<T>& R) {
    incre_nonviscous_force = R;
    trial_nonviscous_force = current_nonviscous_force + incre_nonviscous_force;
}

template<sp_d T> void Factory<T>::update_incre_inertial_force(const Col<T>& R) {
    incre_inertial_force = R;
    trial_inertial_force = current_inertial_force + incre_inertial_force;
}

template<sp_d T> void Factory<T>::update_incre_displacement(const Col<T>& D) {
    incre_displacement = D;
    trial_displacement = current_displacement + incre_displacement;
}

template<sp_d T> void Factory<T>::update_incre_velocity(const Col<T>& V) {
    incre_velocity = V;
    trial_velocity = current_velocity + incre_velocity;
}

template<sp_d T> void Factory<T>::update_incre_acceleration(const Col<T>& A) {
    incre_acceleration = A;
    trial_acceleration = current_acceleration + incre_acceleration;
}

template<sp_d T> void Factory<T>::update_incre_temperature(const Col<T>& M) {
    incre_temperature = M;
    trial_temperature = current_temperature + incre_temperature;
}

template<sp_d T> void Factory<T>::update_current_time(const T M) {
    trial_time = current_time = M;
    incre_time = T(0);
}

template<sp_d T> void Factory<T>::update_current_load_factor(const Col<T>& L) {
    trial_load_factor = current_load_factor = L;
    incre_load_factor.zeros();
}

template<sp_d T> void Factory<T>::update_current_load(const Col<T>& L) {
    trial_load = current_load = L;
    incre_load.zeros();
}

template<sp_d T> void Factory<T>::update_current_settlement(const Col<T>& S) {
    trial_settlement = current_settlement = S;
    incre_settlement.zeros();
}

template<sp_d T> void Factory<T>::update_current_resistance(const Col<T>& R) {
    trial_resistance = current_resistance = R;
    incre_resistance.zeros();
}

template<sp_d T> void Factory<T>::update_current_damping_force(const Col<T>& R) {
    trial_damping_force = current_damping_force = R;
    incre_damping_force.zeros();
}

template<sp_d T> void Factory<T>::update_current_nonviscous_force(const Col<T>& R) {
    trial_nonviscous_force = current_nonviscous_force = R;
    incre_nonviscous_force.zeros();
}

template<sp_d T> void Factory<T>::update_current_inertial_force(const Col<T>& R) {
    trial_inertial_force = current_inertial_force = R;
    incre_inertial_force.zeros();
}

template<sp_d T> void Factory<T>::update_current_displacement(const Col<T>& D) {
    trial_displacement = current_displacement = D;
    incre_displacement.zeros();
}

template<sp_d T> void Factory<T>::update_current_velocity(const Col<T>& V) {
    trial_velocity = current_velocity = V;
    incre_velocity.zeros();
}

template<sp_d T> void Factory<T>::update_current_acceleration(const Col<T>& A) {
    trial_acceleration = current_acceleration = A;
    incre_acceleration.zeros();
}

template<sp_d T> void Factory<T>::update_current_temperature(const Col<T>& M) {
    trial_temperature = current_temperature = M;
    incre_temperature.zeros();
}

template<sp_d T> void Factory<T>::update_trial_time_by(const T M) {
    trial_time += M;
    incre_time = trial_time - current_time;
}

template<sp_d T> void Factory<T>::update_trial_load_factor_by(const Col<T>& L) {
    trial_load_factor += L;
    incre_load_factor = trial_load_factor - current_load_factor;
}

template<sp_d T> void Factory<T>::update_trial_load_by(const Col<T>& L) {
    trial_load += L;
    incre_load = trial_load - current_load;
}

template<sp_d T> void Factory<T>::update_trial_settlement_by(const Col<T>& S) {
    trial_settlement += S;
    incre_settlement = trial_settlement - current_settlement;
}

template<sp_d T> void Factory<T>::update_trial_resistance_by(const Col<T>& R) {
    trial_resistance += R;
    incre_resistance = trial_resistance - current_resistance;
}

template<sp_d T> void Factory<T>::update_trial_damping_force_by(const Col<T>& R) {
    trial_damping_force += R;
    incre_damping_force = trial_damping_force - current_damping_force;
}

template<sp_d T> void Factory<T>::update_trial_nonviscous_force_by(const Col<T>& R) {
    trial_nonviscous_force += R;
    incre_nonviscous_force = trial_nonviscous_force - current_nonviscous_force;
}

template<sp_d T> void Factory<T>::update_trial_inertial_force_by(const Col<T>& R) {
    trial_inertial_force += R;
    incre_inertial_force = trial_inertial_force - current_inertial_force;
}

template<sp_d T> void Factory<T>::update_trial_displacement_by(const Col<T>& D) {
    trial_displacement += D;
    incre_displacement = trial_displacement - current_displacement;
}

template<sp_d T> void Factory<T>::update_trial_velocity_by(const Col<T>& V) {
    trial_velocity += V;
    incre_velocity = trial_velocity - current_velocity;
}

template<sp_d T> void Factory<T>::update_trial_acceleration_by(const Col<T>& A) {
    trial_acceleration += A;
    incre_acceleration = trial_acceleration - current_acceleration;
}

template<sp_d T> void Factory<T>::update_trial_temperature_by(const Col<T>& M) {
    trial_temperature += M;
    incre_temperature = trial_temperature - current_temperature;
}

template<sp_d T> void Factory<T>::update_incre_time_by(const T M) {
    incre_time += M;
    trial_time = current_time + incre_time;
}

template<sp_d T> void Factory<T>::update_incre_load_factor_by(const Col<T>& L) {
    incre_load_factor += L;
    trial_load_factor = current_load_factor + incre_load_factor;
}

template<sp_d T> void Factory<T>::update_incre_load_by(const Col<T>& L) {
    incre_load += L;
    trial_load = current_load + incre_load;
}

template<sp_d T> void Factory<T>::update_incre_settlement_by(const Col<T>& S) {
    incre_settlement += S;
    trial_settlement = current_settlement + incre_settlement;
}

template<sp_d T> void Factory<T>::update_incre_resistance_by(const Col<T>& R) {
    incre_resistance += R;
    trial_resistance = current_resistance + incre_resistance;
}

template<sp_d T> void Factory<T>::update_incre_damping_force_by(const Col<T>& R) {
    incre_damping_force += R;
    trial_damping_force = current_damping_force + incre_damping_force;
}

template<sp_d T> void Factory<T>::update_incre_nonviscous_force_by(const Col<T>& R) {
    incre_nonviscous_force += R;
    trial_nonviscous_force = current_nonviscous_force + incre_nonviscous_force;
}

template<sp_d T> void Factory<T>::update_incre_inertial_force_by(const Col<T>& R) {
    incre_inertial_force += R;
    trial_inertial_force = current_inertial_force + incre_inertial_force;
}

template<sp_d T> void Factory<T>::update_incre_displacement_by(const Col<T>& D) {
    incre_displacement += D;
    trial_displacement = current_displacement + incre_displacement;
}

template<sp_d T> void Factory<T>::update_incre_velocity_by(const Col<T>& V) {
    incre_velocity += V;
    trial_velocity = current_velocity + incre_velocity;
}

template<sp_d T> void Factory<T>::update_incre_acceleration_by(const Col<T>& A) {
    incre_acceleration += A;
    trial_acceleration = current_acceleration + incre_acceleration;
}

template<sp_d T> void Factory<T>::update_incre_temperature_by(const Col<T>& M) {
    incre_temperature += M;
    trial_temperature = current_temperature + incre_temperature;
}

template<sp_d T> void Factory<T>::update_current_time_by(const T M) {
    trial_time = current_time += M;
    incre_time = 0.;
}

template<sp_d T> void Factory<T>::update_current_load_factor_by(const Col<T>& L) {
    trial_load_factor = current_load_factor += L;
    incre_load_factor.zeros();
}

template<sp_d T> void Factory<T>::update_current_load_by(const Col<T>& L) {
    trial_load = current_load += L;
    incre_load.zeros();
}

template<sp_d T> void Factory<T>::update_current_settlement_by(const Col<T>& S) {
    trial_settlement = current_settlement += S;
    incre_settlement.zeros();
}

template<sp_d T> void Factory<T>::update_current_resistance_by(const Col<T>& R) {
    trial_resistance = current_resistance += R;
    incre_resistance.zeros();
}

template<sp_d T> void Factory<T>::update_current_damping_force_by(const Col<T>& R) {
    trial_damping_force = current_damping_force += R;
    incre_damping_force.zeros();
}

template<sp_d T> void Factory<T>::update_current_nonviscous_force_by(const Col<T>& R) {
    trial_nonviscous_force = current_nonviscous_force += R;
    incre_nonviscous_force.zeros();
}

template<sp_d T> void Factory<T>::update_current_inertial_force_by(const Col<T>& R) {
    trial_inertial_force = current_inertial_force += R;
    incre_inertial_force.zeros();
}

template<sp_d T> void Factory<T>::update_current_displacement_by(const Col<T>& D) {
    trial_displacement = current_displacement += D;
    incre_displacement.zeros();
}

template<sp_d T> void Factory<T>::update_current_velocity_by(const Col<T>& V) {
    trial_velocity = current_velocity += V;
    incre_velocity.zeros();
}

template<sp_d T> void Factory<T>::update_current_acceleration_by(const Col<T>& A) {
    trial_acceleration = current_acceleration += A;
    incre_acceleration.zeros();
}

template<sp_d T> void Factory<T>::update_current_temperature_by(const Col<T>& M) {
    trial_temperature = current_temperature += M;
    incre_temperature.zeros();
}

#endif // FACTORY_HPP

//! @}
