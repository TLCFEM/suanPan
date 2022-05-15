/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class Frequency
 * @brief A Frequency class.
 * @author tlc
 * @date 06/07/2018
 * @version 0.2.0
 * @file Frequency.h
 * @addtogroup Step
 * @{
 */

#ifndef FREQUENCY_H
#define FREQUENCY_H

#include <Step/Step.h>

class Frequency final : public Step {
    const unsigned eigen_number;

public:
    explicit Frequency(unsigned = 0, unsigned = 4);

    int initialize() override;

    int analyze() override;

    void set_eigen_number(unsigned) const;
    [[nodiscard]] unsigned get_eigen_number() const;
};

#endif

//! @}
