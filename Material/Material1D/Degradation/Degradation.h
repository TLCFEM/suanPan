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
 * @class Degradation
 * @brief The Degradation class.
 *
 * @author tlc
 * @date 11/05/2019
 * @version 0.1.0
 * @file Degradation.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef DEGRADATION_H
#define DEGRADATION_H

#include <Material/Material1D/Material1D.h>

class Degradation : public Material1D {
    const unsigned mat_tag;

protected:
    unique_ptr<Material> base;

    [[nodiscard]] virtual vec compute_positive_degradation(double) const = 0; // positive region
    [[nodiscard]] virtual vec compute_negative_degradation(double) const = 0; // negative region

public:
    Degradation(
        unsigned, // tag
        unsigned  // material tag
    );
    Degradation(const Degradation&);
    Degradation(Degradation&&) noexcept = delete;
    Degradation& operator=(const Degradation&) = delete;
    Degradation& operator=(Degradation&&) noexcept = delete;
    ~Degradation() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

class StrainDegradation : public Degradation {
public:
    using Degradation::Degradation;

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    std::vector<vec> record(OutputType) override;
};

class StressDegradation : public Degradation {
public:
    using Degradation::Degradation;

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    std::vector<vec> record(OutputType) override;
};

#endif

//! @}
