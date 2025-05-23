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
 * @class Subloading1D
 * @brief A Subloading1D material class.
 * @author tlc
 * @date 24/09/2024
 * @version 0.1.0
 * @file Subloading1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef SUBLOADING1D_H
#define SUBLOADING1D_H

#include <Material/Material1D/Material1D.h>

struct DataSubloading1D {
    class Saturation {
        static const double root_one_half;

        const double rate, bound;

    public:
        Saturation(const double R, const double B)
            : rate(R)
            , bound(B) {}

        [[nodiscard]] double r() const { return rate * root_one_half; }

        [[nodiscard]] double b() const { return bound; }

        [[nodiscard]] double rb() const { return r() * b(); }
    };

    const double elastic; // elastic modulus
    const double initial_iso, k_iso, saturation_iso, m_iso;
    const double initial_kin, k_kin, saturation_kin, m_kin;
    const double u;
    const double cv;
    const double mu;
    const double nv;

    const std::vector<Saturation> b, c;
};

class Subloading1D final : protected DataSubloading1D, public Material1D {
    static constexpr unsigned max_iteration = 20u;
    static constexpr double z_bound = 1E-15;
    static const double rate_bound;

    static pod2 yield_ratio(double);

    const double* incre_time = nullptr;

    const bool is_viscous = mu > 0. && nv > 0.;

public:
    Subloading1D(
        unsigned,           // tag
        DataSubloading1D&&, // data
        double = 0.         // density
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
