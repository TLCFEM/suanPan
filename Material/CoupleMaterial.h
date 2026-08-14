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
 * @class CoupleMaterial
 * @brief A CoupleMaterial abstract base class.
 * @author tlc
 * @date 14/08/2026
 * @version 0.1.0
 * @file CoupleMaterial.h
 * @addtogroup Material
 * @{
 */

#ifndef COUPLEMATERIAL_H
#define COUPLEMATERIAL_H

#include <suanPan.h>

class DomainBase;

struct DataCoupleMaterial {
    const double characteristic_length;

    vec current_curvature{};
    vec current_moment{};

    vec trial_curvature{};
    vec trial_moment{};

    vec incre_curvature{};
    vec incre_moment{};

    mat initial_stiffness{};
    mat current_stiffness{};
    mat trial_stiffness{};
};

class CoupleMaterial : protected DataCoupleMaterial {
    friend void ConstantStiffness(CoupleMaterial*);

public:
    explicit CoupleMaterial(double);
    virtual ~CoupleMaterial() = default;

    [[nodiscard]] virtual std::unique_ptr<CoupleMaterial> unique_copy() const = 0;

    virtual void initialize(const shared_ptr<DomainBase>&) {}

    const vec& get_trial_curvature();
    const vec& get_trial_moment();
    const mat& get_trial_stiffness();

    const vec& get_current_curvature();
    const vec& get_current_moment();
    const mat& get_current_stiffness();

    [[nodiscard]] const mat& get_initial_stiffness() const;

    int update_incre_status(double);
    int update_incre_status(double, double);
    int update_incre_status(double, double, double);
    int update_trial_status(double);
    int update_trial_status(double, double);
    int update_trial_status(double, double, double);

    int update_incre_status(const vec&);
    int update_incre_status(const vec&, const vec&);
    int update_incre_status(const vec&, const vec&, const vec&);
    int update_trial_status(const vec&);
    int update_trial_status(const vec&, const vec&);
    int update_trial_status(const vec&, const vec&, const vec&);

    int clear_status();
    int commit_status();
    int reset_status();
};

class IsotropicCouple final : public CoupleMaterial {
    const double elastic_modulus, poissons_ratio;

public:
    IsotropicCouple(double, double, double);

    [[nodiscard]] std::unique_ptr<CoupleMaterial> unique_copy() const override { return std::make_unique<IsotropicCouple>(*this); }

    void initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
