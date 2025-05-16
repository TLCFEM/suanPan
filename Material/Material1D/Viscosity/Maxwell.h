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
 * @class Maxwell
 * @brief A 1D Maxwell material class.
 * @author tlc
 * @date 07/07/2018
 * @version 0.3.0
 * @file Maxwell.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef MAXWELL_H
#define MAXWELL_H

#include <Material/Material1D/Material1D.h>
#include <Toolbox/ResourceHolder.h>

class Maxwell final : public Material1D {
    static constexpr unsigned max_iteration = 20u;

    const double* incre_time = nullptr;

    unsigned counter = 0, delay_counter = 0;
    const unsigned damper_tag, spring_tag;
    const unsigned proceed;

    const bool use_matrix;

    double beta;

    ResourceHolder<Material> damper, spring;

public:
    Maxwell(
        unsigned,     // tag
        unsigned,     // damper tag
        unsigned,     // spring tag
        bool = false, // if to use matrix
        unsigned = 0, // if to process when fails to converge
        double = 0.   // beta
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

//! @}
