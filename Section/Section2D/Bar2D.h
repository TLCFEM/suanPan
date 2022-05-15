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
 * @class Bar2D
 * @brief A Bar2D class.
 * @author tlc
 * @date 05/06/2018
 * @version 0.1.0
 * @file Bar2D.h
 * @addtogroup Section
 * @{
 */

#ifndef BAR2D_H
#define BAR2D_H

#include <Section/Section2D/Section2D.h>

class Bar2D final : public Section2D {
    unique_ptr<Material> s_material;

public:
    Bar2D(unsigned,   // tag
          double,     // area
          unsigned,   // material tag
          double = 0. // eccentricity
    );
    Bar2D(const Bar2D&);
    Bar2D(Bar2D&&) noexcept = delete;            // move forbidden
    Bar2D& operator=(const Bar2D&) = delete;     // assign forbidden
    Bar2D& operator=(Bar2D&&) noexcept = delete; // assign forbidden
    ~Bar2D() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
