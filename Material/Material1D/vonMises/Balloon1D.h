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
 * @class Balloon1D
 * @brief A Balloon1D material class.
 * @author tlc
 * @date 27/09/2025
 * @version 0.1.0
 * @file Balloon1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BALLOON1D_H
#define BALLOON1D_H

#include <Material/Material1D/Material1D.h>

struct DataBalloon1D {
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
    const double initial_iso_m, k_iso_m, saturation_iso_m, m_iso_m;
    const double initial_iso_r, k_iso_r, saturation_iso_r, m_iso_r;
    const double initial_kin, k_kin, saturation_kin, m_kin;
    const double u;

    const unsigned zr_size;

    const std::vector<Saturation> b, c;
};

class Balloon1D final : protected DataBalloon1D, public Material1D {
    static constexpr unsigned max_iteration = 20u;
    static constexpr double z_bound = 1E-15;
    static const double rate_bound;

    static pod2 yield_ratio(double);

    class memory {
        std::vector<double> buffer;
        std::size_t head;

    public:
        memory(const std::size_t size)
            : buffer(size, 0.)
            , head(0) {}

        auto& enqueue(const double value) {
            buffer[head] = value;
            head = (head + 1) % buffer.size();
            return *this;
        }

        auto zeros() {
            std::fill(buffer.begin(), buffer.end(), 0.);
            head = 0;
        }

        auto max() { return *std::max_element(buffer.cbegin(), buffer.cend()); }
        auto min() { return *std::min_element(buffer.cbegin(), buffer.cend()); }
        auto mean() { return std::accumulate(buffer.cbegin(), buffer.cend(), 0.) / buffer.size(); }
    };

    memory current_zr{zr_size}, trial_zr{zr_size};

    auto push_zr(const double zr) { return (trial_zr = current_zr).enqueue(zr).mean(); }

public:
    Balloon1D(
        unsigned,        // tag
        DataBalloon1D&&, // data
        double = 0.      // density
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
