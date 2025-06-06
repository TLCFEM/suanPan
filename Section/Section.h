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
 * @class Section
 * @brief A Section class.
 * @author tlc
 * @date 06/07/2018
 * @version 0.2.0
 * @file Section.h
 * @addtogroup Section
 * @{
 */

#ifndef SECTION_H
#define SECTION_H

#include <Domain/Tag.h>

enum class SectionType : unsigned {
    D0 = 0,
    D1 = 1,
    D2 = 2,
    D3 = 3,
    NM2D = 3,
    NM3D = 6,
    OS3D = 12,
};

enum class OutputType;

class DomainBase;
class Material;

struct DataSection {
    const unsigned material_tag; // material tag

    const SectionType section_type;

    const vec eccentricity;

    const double area;
    const double linear_density = 0.;
    const double characteristic_length = -1.;

    vec trial_deformation{};   // trial deformation
    vec current_deformation{}; // current deformation

    vec trial_deformation_rate{};   // trial deformation rate
    vec current_deformation_rate{}; // current deformation rate

    vec trial_resistance{};   // trial resistance
    vec current_resistance{}; // current resistance

    mat initial_stiffness{}; // stiffness matrix
    mat current_stiffness{}; // stiffness matrix
    mat trial_stiffness{};   // stiffness matrix

    mat initial_geometry{}; // geometry matrix
    mat current_geometry{}; // geometry matrix
    mat trial_geometry{};   // geometry matrix
};

class Section : protected DataSection, public CopiableTag {
    const bool initialized = false;
    const bool symmetric = false;

public:
    explicit Section(
        unsigned = 0,                  // section tag
        SectionType = SectionType::D0, // section type
        unsigned = 0,                  // material tag
        double = 0.,                   // area
        vec&& = {0., 0.}               // eccentricity
    );

    [[nodiscard]] SectionType get_section_type() const;
    [[nodiscard]] double get_area() const;
    [[nodiscard]] double get_linear_density() const;

    int initialize_base(const shared_ptr<DomainBase>&);

    virtual int initialize(const shared_ptr<DomainBase>&) = 0;

    void set_initialized(bool) const;
    void set_symmetric(bool) const;
    [[nodiscard]] bool is_initialized() const;
    [[nodiscard]] bool is_symmetric() const;

    void set_eccentricity(const vec&) const;
    [[nodiscard]] const vec& get_eccentricity() const;

    virtual void set_characteristic_length(double) const;
    [[nodiscard]] double get_characteristic_length() const;

    [[nodiscard]] virtual const vec& get_trial_deformation() const;
    [[nodiscard]] virtual const vec& get_trial_deformation_rate() const;
    [[nodiscard]] virtual const vec& get_trial_resistance() const;
    [[nodiscard]] virtual const mat& get_trial_stiffness() const;
    [[nodiscard]] virtual const mat& get_trial_geometry() const;

    [[nodiscard]] virtual const vec& get_current_deformation() const;
    [[nodiscard]] virtual const vec& get_current_deformation_rate() const;
    [[nodiscard]] virtual const vec& get_current_resistance() const;
    [[nodiscard]] virtual const mat& get_current_stiffness() const;
    [[nodiscard]] virtual const mat& get_current_geometry() const;

    [[nodiscard]] virtual const mat& get_initial_stiffness() const;
    [[nodiscard]] virtual const mat& get_initial_geometry() const;

    virtual unique_ptr<Section> get_copy() = 0;

    int update_incre_status(double);
    int update_incre_status(double, double);
    int update_trial_status(double);
    int update_trial_status(double, double);

    virtual int update_incre_status(const vec&);
    virtual int update_incre_status(const vec&, const vec&);
    virtual int update_trial_status(const vec&);
    virtual int update_trial_status(const vec&, const vec&);

    virtual int clear_status() = 0;
    virtual int commit_status() = 0;
    virtual int reset_status() = 0;

    virtual std::vector<vec> record(OutputType);
};

namespace suanpan {
    unique_ptr<Section> make_copy(const shared_ptr<Section>&);
} // namespace suanpan

#endif

//! @}
