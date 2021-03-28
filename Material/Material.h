/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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
 * @class Material
 * @brief A Material abstract base class.
 * @author tlc
 * @date 30/05/2020
 * @version 0.1.2
 * @file Material.h
 * @addtogroup Material
 * @{
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include <Domain/Tag.h>
#include <Section/ParameterType.h>

enum class MaterialType : unsigned { D0 = 0, D1 = 1, D2 = 3, D3 = 6 };

class DomainBase;
enum class OutputType;

using std::vector;

struct MaterialData {
	const double tolerance = 1E-14;
	const double density = 0.; // density
	const MaterialType material_type = MaterialType::D0;

	vec current_strain;      // current status
	vec current_strain_rate; // current status
	vec current_strain_acc;  // current status
	vec current_stress;      // current status
	// vec current_stress_rate; // current status

	vec trial_strain;      // trial status
	vec trial_strain_rate; // trial status
	vec trial_strain_acc;  // trial status
	vec trial_stress;      // trial status
	// vec trial_stress_rate; // trial status

	vec incre_strain;      // incremental status
	vec incre_strain_rate; // incremental status
	vec incre_strain_acc;  // incremental status
	vec incre_stress;      // incremental status
	// vec incre_stress_rate; // incremental status

	vec initial_history; // initial status
	vec current_history; // current status
	vec trial_history;   // trial status

	mat initial_stiffness; // stiffness matrix
	mat current_stiffness; // stiffness matrix
	mat trial_stiffness;   // stiffness matrix

	mat initial_damping; // damping matrix
	mat current_damping; // damping matrix
	mat trial_damping;   // damping matrix

	mat initial_inertial; // inertial matrix
	mat current_inertial; // inertial matrix
	mat trial_inertial;   // inertial matrix
};

class Material : protected MaterialData, public Tag {
	const bool initialized = false;
	const bool symmetric = false;

	friend void ConstantStiffness(MaterialData*);
	friend void ConstantDamping(MaterialData*);
	friend void ConstantInertial(MaterialData*);
	friend void PureWrapper(MaterialData*);
public:
	explicit Material(unsigned = 0,                    // tag
	                  MaterialType = MaterialType::D0, // material type
	                  double = 0.                      // density
	);
	Material(const Material&) = default;
	Material(Material&&) = delete;                 // move forbidden
	Material& operator=(const Material&) = delete; // assign forbidden
	Material& operator=(Material&&) = delete;      // assign forbidden

	~Material() override;

	virtual void initialize(const shared_ptr<DomainBase>&) = 0;

	virtual void initialize_history(unsigned);
	virtual void set_initial_history(const vec&);

	void set_initialized(bool) const;
	void set_symmetric(bool) const;
	[[nodiscard]] bool is_initialized() const;
	[[nodiscard]] bool is_symmetric() const;

	[[nodiscard]] MaterialType get_material_type() const;

	[[nodiscard]] virtual double get_parameter(ParameterType) const;

	virtual const vec& get_trial_strain();
	virtual const vec& get_trial_strain_rate();
	virtual const vec& get_trial_strain_acc();
	virtual const vec& get_trial_stress();
	virtual const mat& get_trial_stiffness();
	virtual const mat& get_trial_secant();
	virtual const mat& get_trial_damping();
	virtual const mat& get_trial_inertial();

	virtual const vec& get_current_strain();
	virtual const vec& get_current_strain_rate();
	virtual const vec& get_current_strain_acc();
	virtual const vec& get_current_stress();
	virtual const mat& get_current_stiffness();
	virtual const mat& get_current_secant();
	virtual const mat& get_current_damping();
	virtual const mat& get_current_inertial();

	[[nodiscard]] virtual const vec& get_initial_history() const;
	[[nodiscard]] virtual const mat& get_initial_stiffness() const;
	[[nodiscard]] virtual const mat& get_initial_damping() const;
	[[nodiscard]] virtual const mat& get_initial_inertial() const;

	virtual unique_ptr<Material> get_copy() = 0;

	int update_incre_status(double);
	int update_incre_status(double, double);
	int update_incre_status(double, double, double);
	int update_trial_status(double);
	int update_trial_status(double, double);
	int update_trial_status(double, double, double);

	virtual int update_incre_status(const vec&);
	virtual int update_incre_status(const vec&, const vec&);
	virtual int update_incre_status(const vec&, const vec&, const vec&);
	virtual int update_trial_status(const vec&);
	virtual int update_trial_status(const vec&, const vec&);
	virtual int update_trial_status(const vec&, const vec&, const vec&);

	virtual int clear_status() = 0;
	virtual int commit_status() = 0;
	virtual int reset_status() = 0;

	virtual vector<vec> record(OutputType);
};

namespace suanpan {
	unique_ptr<Material> make_copy(const shared_ptr<Material>&);
	unique_ptr<Material> make_copy(const unique_ptr<Material>&);
}

#endif

//! @}
