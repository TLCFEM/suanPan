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
 * @class MaterialBase
 * @brief A MaterialBase abstract base class.
 * @author tlc
 * @date 29/05/2020
 * @version 0.1.0
 * @file MaterialBase.h
 * @addtogroup Material
 * @{
 */

#ifndef MATERIALBASE_H
#define MATERIALBASE_H

#include <Domain/Tag.h>

class DomainBase;
enum class MaterialType;
enum class OutputType;
enum class ParameterType;

using std::vector;

class MaterialBase : public Tag {
public:
	explicit MaterialBase(const unsigned T = 0)
		: Tag(T) {}

	MaterialBase(const MaterialBase&) = default;           // allow copy
	MaterialBase(MaterialBase&&) = delete;                 // move forbidden
	MaterialBase& operator=(const MaterialBase&) = delete; // assign forbidden
	MaterialBase& operator=(MaterialBase&&) = delete;      // assign forbidden

	~MaterialBase() override = default;

	virtual void initialize(const shared_ptr<DomainBase>&) = 0;

	virtual void initialize_history(unsigned) = 0;
	virtual void set_initial_history(const vec&) = 0;

	virtual void set_initialized(bool) const = 0;
	virtual void set_symmetric(bool) const = 0;
	[[nodiscard]] virtual bool is_initialized() const = 0;
	[[nodiscard]] virtual bool is_symmetric() const = 0;

	[[nodiscard]] virtual MaterialType get_material_type() const = 0;

	[[nodiscard]] virtual double get_parameter(ParameterType) const = 0;

	[[nodiscard]] virtual const vec& get_trial_strain() = 0;
	[[nodiscard]] virtual const vec& get_trial_strain_rate() = 0;
	[[nodiscard]] virtual const vec& get_trial_strain_acc() = 0;
	[[nodiscard]] virtual const vec& get_trial_stress() = 0;
	[[nodiscard]] virtual const mat& get_trial_stiffness() = 0;
	[[nodiscard]] virtual const mat& get_trial_secant() = 0;
	[[nodiscard]] virtual const mat& get_trial_damping() = 0;
	[[nodiscard]] virtual const mat& get_trial_inertial() = 0;

	[[nodiscard]] virtual const vec& get_current_strain() = 0;
	[[nodiscard]] virtual const vec& get_current_strain_rate() = 0;
	[[nodiscard]] virtual const vec& get_current_strain_acc() = 0;
	[[nodiscard]] virtual const vec& get_current_stress() = 0;
	[[nodiscard]] virtual const mat& get_current_stiffness() = 0;
	[[nodiscard]] virtual const mat& get_current_secant() = 0;
	[[nodiscard]] virtual const mat& get_current_damping() = 0;
	[[nodiscard]] virtual const mat& get_current_inertial() = 0;

	[[nodiscard]] virtual const vec& get_initial_history() const = 0;
	[[nodiscard]] virtual const mat& get_initial_stiffness() const = 0;
	[[nodiscard]] virtual const mat& get_initial_damping() const = 0;
	[[nodiscard]] virtual const mat& get_initial_inertial() const = 0;

	virtual unique_ptr<MaterialBase> get_copy() = 0;

	virtual int update_incre_status(double) = 0;
	virtual int update_incre_status(double, double) = 0;
	virtual int update_incre_status(double, double, double) = 0;
	virtual int update_trial_status(double) = 0;
	virtual int update_trial_status(double, double) = 0;
	virtual int update_trial_status(double, double, double) = 0;

	virtual int update_incre_status(const vec&) = 0;
	virtual int update_incre_status(const vec&, const vec&) = 0;
	virtual int update_incre_status(const vec&, const vec&, const vec&) = 0;
	virtual int update_trial_status(const vec&) = 0;
	virtual int update_trial_status(const vec&, const vec&) = 0;
	virtual int update_trial_status(const vec&, const vec&, const vec&) = 0;

	virtual int clear_status() = 0;
	virtual int commit_status() = 0;
	virtual int reset_status() = 0;

	virtual vector<vec> record(OutputType) = 0;
};

namespace suanpan {
	unique_ptr<MaterialBase> make_copy(const shared_ptr<MaterialBase>&);
	unique_ptr<MaterialBase> make_copy(const unique_ptr<MaterialBase>&);
} // namespace suanpan

#endif

//! @}
