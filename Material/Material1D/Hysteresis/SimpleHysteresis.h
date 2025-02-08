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
 * @class SimpleHysteresis
 * @brief A SimpleHysteresis material class.
 *
 * The SimpleHysteresis class defines a uniaxial hysteresis model with simple hysteresis behaviour.
 *
 * @author tlc
 * @date 28/01/2019
 * @version 0.1.0
 * @file SimpleHysteresis.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef SIMPLEHYSTERESIS_H
#define SIMPLEHYSTERESIS_H

#include <Material/Material1D/Material1D.h>

struct DataSimpleHysteresis {
    const double middle_point;
};

class SimpleHysteresis : protected DataSimpleHysteresis, public Material1D {
    enum class Status {
        NONE,
        CBACKBONE,
        TBACKBONE,
        CINNER,
        TINNER
    };

    Status trial_flag = Status::NONE, current_flag = Status::NONE;

    /**
     * \brief to compute the initial reverse point on the opposite branch
     * \return the reverse point (strain, stress)
     */
    [[nodiscard]] virtual podarray<double> compute_compression_initial_reverse() const = 0;
    /**
     * \brief to compute the initial reverse point on the opposite branch
     * \return the reverse point (strain, stress)
     */
    [[nodiscard]] virtual podarray<double> compute_tension_initial_reverse() const = 0;

    [[nodiscard]] virtual podarray<double> compute_compression_backbone(double) const = 0;
    [[nodiscard]] virtual podarray<double> compute_tension_backbone(double) const = 0;
    [[nodiscard]] virtual double compute_compression_residual(double, double) const = 0;
    [[nodiscard]] virtual double compute_tension_residual(double, double) const = 0;
    [[nodiscard]] podarray<double> compute_compression_inner(double) const;
    [[nodiscard]] podarray<double> compute_tension_inner(double) const;

public:
    SimpleHysteresis(
        unsigned,   // tag
        double,     // middle point
        double = 0. // density
    );

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
