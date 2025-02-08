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
 * @class CustomAmplitude
 * @brief An Amplitude class that can generate Amplitude pattern.
 *
 * @author tlc
 * @date 26/03/2023
 * @version 0.1.0
 * @file CustomAmplitude.h
 * @addtogroup Amplitude
 * @{
 */

#ifndef CUSTOMAMPLITUDE_H
#define CUSTOMAMPLITUDE_H

#include <Load/Amplitude/Amplitude.h>
#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomAmplitude final : public Amplitude {
    const unsigned e_tag;

    ResourceHolder<Expression> expression;

public:
    CustomAmplitude(unsigned, unsigned, unsigned);

    void initialize(const shared_ptr<DomainBase>&) override;

    double get_amplitude(double) override;

    void print() override;
};

#endif

//! @}
