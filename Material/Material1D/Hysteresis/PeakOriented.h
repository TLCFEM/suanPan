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
 * @class PeakOriented
 * @brief A PeakOriented material class.
 * @author tlc
 * @date 13/04/2020
 * @version 0.1.0
 * @file PeakOriented.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef PEAKORIENTED_H
#define PEAKORIENTED_H

#include "SimpleHysteresis.h"

class PeakOriented : public SimpleHysteresis {
    [[nodiscard]] double compute_compression_residual(double, double) const override;
    [[nodiscard]] double compute_tension_residual(double, double) const override;

public:
    PeakOriented(int, double);

    void print() override;
};

#endif

//! @}
