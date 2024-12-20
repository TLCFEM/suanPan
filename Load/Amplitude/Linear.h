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
 * @class Linear
 * @brief A Linear class that can generate Amplitude pattern.
 *
 * @author tlc
 * @date 26/07/2018
 * @version 0.1.0
 * @file Linear.h
 * @addtogroup Amplitude
 * @{
 */

#ifndef LINEAR_H
#define LINEAR_H

#include <Load/Amplitude/Amplitude.h>

class Linear final : public Amplitude {
    const double slope;

public:
    explicit Linear(unsigned, double, unsigned);

    double get_amplitude(double) override;

    void print() override;
};

#endif

//! @}
