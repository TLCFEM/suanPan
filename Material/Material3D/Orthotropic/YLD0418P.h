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
 * @class YLD0418P
 * @brief The YLD0418P class.
 *
 * @author tlc
 * @date 20/07/2025
 * @version 1.0.0
 * @file YLD0418P.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef YLD0418P_H
#define YLD0418P_H

#include <Material/Material3D/Material3D.h>
#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

struct DataYLD0418P {
    const vec modulus, ratio, parameter;
    const double exponent;
    const double ref_stress; // reference stress for scaling
};

class YLD0418P final : protected DataYLD0418P, public Material3D {
    static constexpr double two_third = 2. / 3.;
    static const double root_two_third;
    static constexpr unsigned max_iteration = 20u;
    static constexpr uword sa{0};
    static const span sb;

    const unsigned hardening_tag;

    ResourceHolder<Expression> hardening_expression;

    mat66 dev_ini_stiffness;
    mat66 C1, C2;

    [[nodiscard]] std::pair<double, double> compute_hardening(const double pe) const {
        if(const auto k = hardening_expression->evaluate(pe).at(0); k >= 0.) return {k, hardening_expression->gradient(pe).at(0)};

        return {0., 0.};
    };

    using yield_t = std::tuple<double, vec6, mat66>;

    [[nodiscard]] yield_t compute_yield_surface(const vec3&, const mat33&, const vec3&, const mat33&) const;

public:
    YLD0418P(
        unsigned, // tag
        vec&&,    // elastic modulus
        vec&&,    // poissons ratio
        vec&&,    // parameter (18)
        double,   // exponent
        double,   // reference stress
        unsigned, // hardening expression tag
        double    // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
