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
 * @class SectionShell
 * @brief A SectionShell class.
 * @author tlc
 * @date 18/10/2019
 * @version 0.1.0
 * @file SectionShell.h
 * @addtogroup SectionShell
 * @{
 */

#ifndef SECTIONSHELL_H
#define SECTIONSHELL_H

#include <Domain/Tag.h>
#include <Section/ParameterType.h>

enum class OutputType;

class DomainBase;
class Material;

using std::vector;

struct SectionShellData {
    const unsigned material_tag; // material tag

    const vec eccentricity;

    vec trial_membrane_strain{};
    vec current_membrane_strain{};
    vec trial_plate_strain{};
    vec current_plate_strain{};

    vec trial_membrane_strain_rate{};
    vec current_membrane_strain_rate{};
    vec trial_plate_strain_rate{};
    vec current_plate_strain_rate{};

    vec trial_membrane_stress{};
    vec current_membrane_stress{};
    vec trial_plate_stress{};
    vec current_plate_stress{};

    mat initial_membrane_stiffness{};
    mat trial_membrane_stiffness{};
    mat current_membrane_stiffness{};
    mat initial_plate_stiffness{};
    mat current_plate_stiffness{};
    mat trial_plate_stiffness{};
};

class SectionShell : protected SectionShellData, public Tag {
    const bool symmetric = false;
    const bool initialized = false;

public:
    explicit SectionShell(unsigned = 0,    // section tag
                          unsigned = 0,    // material tag
                          vec&& = {0., 0.} // eccentricity
    );
    SectionShell(const SectionShell&) = default;           // default copy ctor
    SectionShell(SectionShell&&) = delete;                 // move forbidden
    SectionShell& operator=(const SectionShell&) = delete; // assign forbidden
    SectionShell& operator=(SectionShell&&) = delete;      // assign forbidden

    ~SectionShell() override = default;

    virtual int initialize(const shared_ptr<DomainBase>&) = 0;

    void set_initialized(bool) const;
    void set_symmetric(bool) const;
    [[nodiscard]] bool is_initialized() const;
    [[nodiscard]] bool is_symmetric() const;

    void set_eccentricity(const vec&) const;
    [[nodiscard]] const vec& get_eccentricity() const;

    [[nodiscard]] virtual const vec& get_trial_membrane_strain() const;
    [[nodiscard]] virtual const vec& get_trial_membrane_strain_rate() const;
    [[nodiscard]] virtual const vec& get_trial_plate_strain() const;
    [[nodiscard]] virtual const vec& get_trial_plate_strain_rate() const;
    [[nodiscard]] virtual const vec& get_trial_membrane_stress() const;
    [[nodiscard]] virtual const vec& get_trial_plate_stress() const;
    [[nodiscard]] virtual const mat& get_trial_membrane_stiffness() const;
    [[nodiscard]] virtual const mat& get_trial_plate_stiffness() const;

    [[nodiscard]] virtual const vec& get_current_membrane_strain() const;
    [[nodiscard]] virtual const vec& get_current_membrane_strain_rate() const;
    [[nodiscard]] virtual const vec& get_current_plate_strain() const;
    [[nodiscard]] virtual const vec& get_current_plate_strain_rate() const;
    [[nodiscard]] virtual const vec& get_current_membrane_stress() const;
    [[nodiscard]] virtual const vec& get_current_plate_stress() const;
    [[nodiscard]] virtual const mat& get_current_membrane_stiffness() const;
    [[nodiscard]] virtual const mat& get_current_plate_stiffness() const;

    [[nodiscard]] virtual const mat& get_initial_membrane_stiffness() const;
    [[nodiscard]] virtual const mat& get_initial_plate_stiffness() const;

    virtual unique_ptr<SectionShell> get_copy() = 0;

    virtual double get_parameter(ParameterType);

    int update_incre_status(double, double);
    int update_incre_status(double, double, double, double);
    int update_trial_status(double, double);
    int update_trial_status(double, double, double, double);

    virtual int update_incre_status(const vec&, const vec&);
    virtual int update_incre_status(const vec&, const vec&, const vec&, const vec&);
    virtual int update_trial_status(const vec&, const vec&);
    virtual int update_trial_status(const vec&, const vec&, const vec&, const vec&);

    virtual int clear_status() = 0;
    virtual int commit_status() = 0;
    virtual int reset_status() = 0;

    virtual vector<vec> record(OutputType);
};

namespace suanpan {
    unique_ptr<SectionShell> make_copy(const shared_ptr<SectionShell>&);
    unique_ptr<SectionShell> make_copy(const unique_ptr<SectionShell>&);
} // namespace suanpan

#endif

//! @}
