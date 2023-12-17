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
 * @class Static
 * @brief A Static class.
 *
 * This class corresponds to the Static, General step in Abaqus, which handles a static problem using Newton (or quasi Newton) solvers.
 *
 * @author tlc
 * @date 26/05/2018
 * @version 0.1.1
 * @file Static.h
 * @addtogroup Step
 * @{
 */

#ifndef STATIC_H
#define STATIC_H

#include <Step/Step.h>

class Static : public Step {
public:
    explicit Static(
        unsigned = 0, // tag
        double = 1.   // step time period
    );

    int initialize() override;

    int analyze() override;
};

#endif

//! @}
