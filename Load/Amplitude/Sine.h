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
 * @class Sine
 * @brief An Amplitude class that can generate Amplitude pattern.
 *
 * @author tlc
 * @date 26/10/2017
 * @version 0.1.0
 * @file Sine.h
 * @addtogroup Amplitude
 * @{
 */

#ifndef SINE_H
#define SINE_H

#include <Load/Amplitude/Amplitude.h>

class Sine final : public Amplitude {
    const double period;

    const std::vector<double> amp;

public:
    Sine(unsigned, double, std::vector<double>&&);

    unique_ptr<Amplitude> get_copy() override;

    double get_amplitude(double) override;

    void print() override;
};

#endif

//! @}
