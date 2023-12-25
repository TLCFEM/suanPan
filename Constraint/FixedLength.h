/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class FixedLength
 * @brief A FixedLength class.
 *
 * The `FixedLength` constraint applies constraint to two nodes so that the
 * distance remain constant between those two nodes.
 *
 * @author tlc
 * @date 07/03/2021
 * @version 0.1.0
 * @file FixedLength.h
 * @addtogroup Constraint
 * @{
 */

#ifndef FIXEDLENGTH_H
#define FIXEDLENGTH_H

#include "Constraint.h"

class FixedLength : public Constraint {
    vec coor;

protected:
    const bool min_bound = false, max_bound = false;
    const double min_gap = 0., max_gap = 0.;

public:
    FixedLength(unsigned, unsigned, unsigned, uvec&&);

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;

    void update_status(const vec&) override;
    void commit_status() override;
    void clear_status() override;
    void reset_status() override;
};

class MinimumGap final : public FixedLength {
public:
    MinimumGap(unsigned, unsigned, unsigned, double, uvec&&);
};

class MaximumGap final : public FixedLength {
public:
    MaximumGap(unsigned, unsigned, unsigned, double, uvec&&);
};

class Sleeve final : public FixedLength {
public:
    Sleeve(unsigned, unsigned, unsigned, double, double, uvec&&);
};

#endif

//! @}
