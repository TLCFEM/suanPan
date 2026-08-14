/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
#include <array>

using pod2 = std::array<double, 2>;
using pod6 = std::array<double, 6>;

enum class MaterialType : unsigned {
    D0 = 0,
    D1 = 1,
    D2 = 3,
    D3 = 6,
    DS = 10,
    OS = 3
};

enum class PlaneType : unsigned {
    S = 1,
    E = 2,
    A = 3,
    N = 0
};

class DomainBase;
enum class OutputType;

struct DataMaterial {
    const double density = 0.;
    const MaterialType material_type = MaterialType::D0;
    const PlaneType plane_type = PlaneType::N;

    const double tolerance = 1E-14;

    vec current_strain{}; // current status
    vec trial_strain{};   // trial status
    vec incre_strain{};   // incremental status

    vec current_strain_rate{}; // current status
    vec trial_strain_rate{};   // trial status
    vec incre_strain_rate{};   // incremental status

    vec current_strain_acc{}; // current status
    vec trial_strain_acc{};   // trial status
    vec incre_strain_acc{};   // incremental status

    vec current_stress{}; // current status
    vec trial_stress{};   // trial status
    vec incre_stress{};   // incremental status

    // vec current_stress_rate{}; // current status
    // vec trial_stress_rate{};   // trial status
    // vec incre_stress_rate{};   // incremental status

    vec initial_history{}; // initial status
    vec current_history{}; // current status
    vec trial_history{};   // trial status

    mat initial_stiffness{}; // stiffness matrix
    mat current_stiffness{}; // stiffness matrix
    mat trial_stiffness{};   // stiffness matrix

    mat initial_damping{}; // damping matrix
    mat current_damping{}; // damping matrix
    mat trial_damping{};   // damping matrix

    mat initial_inertial{}; // inertial matrix
    mat current_inertial{}; // inertial matrix
    mat trial_inertial{};   // inertial matrix
};

class Material : protected DataMaterial, public CopyableTag {
    const bool initialized = false;
    const bool symmetric = false;

    friend void ConstantStiffness(DataMaterial*);
    friend void ConstantDamping(DataMaterial*);
    friend void ConstantInertial(DataMaterial*);
    friend void PureWrapper(DataMaterial*);

public:
    enum class Parameter {
        ELASTIC,
        POISSON,
        SHEAR,
        BULK,
        PEAKSTRAIN,
        CRACKSTRAIN
    };

    explicit Material(
        unsigned = 0,                    // tag
        MaterialType = MaterialType::D0, // material type
        double = 0.                      // density
    );

    [[nodiscard]] double get_density() const;
    [[nodiscard]] MaterialType get_material_type() const;
    [[nodiscard]] PlaneType get_plane_type() const;

    int initialize_base(const shared_ptr<DomainBase>&);

    virtual int initialize(const shared_ptr<DomainBase>&) = 0;

    virtual void initialize_history(unsigned);
    virtual void set_initial_history(const vec&);

    void set_initialized(bool) const;
    void set_symmetric(bool) const;
    [[nodiscard]] bool is_initialized() const;
    [[nodiscard]] bool is_symmetric() const;

    [[nodiscard]] virtual double get(Parameter) const;

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

    virtual unique_ptr<Material> unique_copy() = 0;

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

    [[nodiscard]] virtual std::vector<vec> record(OutputType) const;

protected:
    class MaterialProperty {
        const double e, v;

    public:
        MaterialProperty(const double E, const double P)
            : e(E)
            , v(P) {}

        double operator()(const Parameter P) const {
            switch(P) {
            case Parameter::ELASTIC:
                return e;
            case Parameter::POISSON:
                return v;
            case Parameter::SHEAR:
                return e / (2. + 2. * v);
            case Parameter::BULK:
                return e / (3. - 6. * v);
            default:
                return 0.;
            }
        }
    };
};

namespace suanpan {
    unique_ptr<Material> unique_copy(const shared_ptr<Material>&);
} // namespace suanpan

#endif

//! @}
