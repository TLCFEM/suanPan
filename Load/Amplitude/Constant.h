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
 * @class Constant
 * @brief An Amplitude class that can generate Amplitude pattern.
 *
 * \f{gather}{a=1\f}
 *
 * @author tlc
 * @date 03/07/2017
 * @version 0.1.0
 * @file Constant.h
 * @addtogroup Amplitude
 * @{
 */

#ifndef CONSTANT_H
#define CONSTANT_H

#include <Load/Amplitude/Amplitude.h>

class Constant final : public Amplitude {
public:
    using Amplitude::Amplitude;

    unique_ptr<Amplitude> get_copy() override;

    double get_amplitude(double) override;

    void print() override;
};

#endif

//! @}
