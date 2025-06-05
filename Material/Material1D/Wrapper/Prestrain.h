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
 * @class Prestrain
 * @brief A Prestrain material class.
 *
 * The Prestrain material class is a container to hold another uniaxial material.
 * It applies a prestrain to the contained material.
 *
 * @author tlc
 * @date 05/06/2025
 * @version 0.1.0
 * @file Prestrain.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef PRESTRAIN_H
#define PRESTRAIN_H

#include <Load/Amplitude/Amplitude.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/ResourceHolder.h>

class Prestrain final : public Material1D {
    const unsigned base_tag, amplitude_tag;

    const double magnitude;

    const double* analysis_time{};

    std::shared_ptr<Amplitude> amplitude;

    ResourceHolder<Material> base;

    double get_prestrain() const;

public:
    Prestrain(
        unsigned, // tag
        unsigned, // base tag
        unsigned, // amplitude tag
        double    // magnitude
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;
    int update_trial_status(const vec&, const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    std::vector<vec> record(OutputType) override;

    void print() override;
};

#endif
