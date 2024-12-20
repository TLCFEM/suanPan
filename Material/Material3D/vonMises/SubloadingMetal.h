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
 * @class SubloadingMetal
 * @brief A SubloadingMetal material class.
 * @author tlc
 * @date 29/09/2024
 * @version 0.1.0
 * @file SubloadingMetal.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef SUBLOADINGMETAL_H
#define SUBLOADINGMETAL_H

#include <Material/Material3D/Material3D.h>

struct DataSubloadingMetal {
    class Saturation {
        const double rate, bound;

    public:
        Saturation(const double R, const double B)
            : rate(R)
            , bound(B) {}

        [[nodiscard]] double r() const { return rate; }

        [[nodiscard]] double b() const { return bound; }

        [[nodiscard]] double rb() const { return r() * b(); }
    };

    const double elastic; // elastic modulus
    const double poissons_ratio;
    const double initial_iso;
    const double k_iso;
    const double saturation_iso;
    const double m_iso;
    const double initial_kin;
    const double k_kin;
    const double saturation_kin;
    const double m_kin;
    const double u;

    const Saturation b;
    const Saturation c;
};

class SubloadingMetal final : protected DataSubloadingMetal, public Material3D {
    static constexpr unsigned max_iteration = 20u;
    static constexpr double two_third = 2. / 3.;
    static const double root_two_third;
    static constexpr double z_bound = 1E-15;
    static const double rate_bound;
    static const mat unit_dev_tensor;

    static vec2 yield_ratio(double);

    const double double_shear = elastic / (1. + poissons_ratio); // double shear modulus

public:
    SubloadingMetal(
        unsigned,              // tag
        DataSubloadingMetal&&, // data
        double = 0.            // density
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
