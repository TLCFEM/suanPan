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
 * @class ArcLength
 * @brief A ArcLength class.
 *
 * This class corresponds to the Static, Riks step in Abaqus, which handles a
 * static problem using arc-length solvers.
 *
 * @author tlc
 * @date 27/09/2017
 * @version 0.1.2
 * @file ArcLength.h
 * @addtogroup Step
 * @{
 */

#ifndef ARCLENGTH_H
#define ARCLENGTH_H

#include <Step/Step.h>

class ArcLength final : public Step {
public:
    explicit ArcLength(unsigned);

    int initialize() override;

    int analyze() override;
};

#endif // ARCLENGTH_H

//! @}
