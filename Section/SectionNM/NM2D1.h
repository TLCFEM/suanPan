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
 * @class NM2D1
 * @brief A NM2D1 class.
 * @author tlc
 * @date 26/11/2021
 * @version 0.1.0
 * @file NM2D1.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef NM2D1_H
#define NM2D1_H

#include <Section/SectionNM/SectionNM2D.h>

class NM2D1 final : public SectionNM2D {
public:
    NM2D1(
        unsigned, // tag
        double, double, double
    );

    unique_ptr<Section> get_copy() override;

    int update_trial_status(const vec&) override;

    void print() override;
};

#endif

//! @}
