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
 * @class TabularSpline
 * @brief A TabularSpline class that can generate Amplitude pattern.
 *
 * @author tlc
 * @date 15/08/2022
 * @version 0.1.0
 * @file TabularSpline.h
 * @addtogroup Amplitude
 * @{
 */

#ifndef TABULARSPLINE_H
#define TABULARSPLINE_H

#include "Tabular.h"

class TabularSpline final : public Tabular {
    vec dt, dy, m;

public:
    using Tabular::Tabular;

    void initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Amplitude> get_copy() override;

    double get_amplitude(double) override;

    void print() override;
};

#endif

//! @}
