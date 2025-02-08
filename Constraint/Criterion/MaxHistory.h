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
 * @class MaxHistory
 * @brief A MaxHistory class.
 *
 * The MaxHistory class tests if the given variable in each element exceeds the given limit, if so, the element is disabled.
 *
 * @author tlc
 * @date 15/09/2020
 * @version 0.1.0
 * @file MaxHistory.h
 * @addtogroup Criterion
 * @{
 */

#ifndef MAXHISTORY_H
#define MAXHISTORY_H

#include <Constraint/Criterion/Criterion.h>

enum class OutputType;

class MaxHistory final : public Criterion {
    const OutputType history_type;
    const double max_history;

public:
    MaxHistory(
        unsigned,   // tag
        unsigned,   // step tag
        OutputType, // history type
        double      // maximum history
    );

    unique_ptr<Criterion> get_copy() override;

    int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
