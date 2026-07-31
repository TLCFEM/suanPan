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
 * @class HistoryCriterion
 * @brief A HistoryCriterion class.
 *
 * The HistoryCriterion class tests if the given variable in each element exceeds the given limit, if so, the element is disabled.
 *
 * @author tlc
 * @date 15/09/2020
 * @version 0.1.0
 * @file MinMaxHistory.h
 * @addtogroup Criterion
 * @{
 */

#ifndef MINMAXHISTORY_H
#define MINMAXHISTORY_H

#include <Constraint/Criterion/Criterion.h>

enum class OutputType;

class HistoryCriterion : public Criterion {
    [[nodiscard]] virtual bool check(double) const = 0;

    const uvec indices;

protected:
    const double limit;

private:
    const OutputType history_type;

public:
    HistoryCriterion(
        unsigned,   // tag
        OutputType, // history type
        double,     // limit
        uvec&&      // indices
    );

    int process(const shared_ptr<DomainBase>&) final;
};

class MaxHistory final : public HistoryCriterion {
    [[nodiscard]] bool check(double) const override;

public:
    using HistoryCriterion::HistoryCriterion;

    unique_ptr<Criterion> unique_copy() override;
};

class MinHistory final : public HistoryCriterion {
    [[nodiscard]] bool check(double) const override;

public:
    using HistoryCriterion::HistoryCriterion;

    unique_ptr<Criterion> unique_copy() override;
};

#endif

//! @}
