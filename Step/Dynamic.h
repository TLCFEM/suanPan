/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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
 * @class Dynamic
 * @brief A Dynamic class.
 * @author tlc
 * @date 03/09/2017
 * @version 0.1.0
 * @file Dynamic.h
 * @addtogroup Step
 * @{
 */

#ifndef DYNAMIC_H
#define DYNAMIC_H

#include <Step/Step.h>

enum class IntegratorType;

class Dynamic final : public Step {
    const IntegratorType analysis_type;

public:
    Dynamic(unsigned, double, IntegratorType);

    int initialize() override;

    int analyze() override;
};

#endif

//! @}
