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
 * @class Section3D
 * @brief A Section3D class.
 * @author tlc
 * @date 27/10/2017
 * @version 0.1.0
 * @file Section3D.h
 * @addtogroup Section-3D
 * @ingroup Section
 * @{
 */

#ifndef SECTION3D_H
#define SECTION3D_H

#include <Section/Section.h>

using std::vector;

class Section3D : public Section {
protected:
    struct IntegrationPoint {
        double coor_y, coor_z, weight;
        unique_ptr<Material> s_material;
        IntegrationPoint(double, double, double, unique_ptr<Material>&&);
        IntegrationPoint(const IntegrationPoint&);
        IntegrationPoint(IntegrationPoint&&) noexcept = default;
        IntegrationPoint& operator=(const IntegrationPoint&) = delete;
        IntegrationPoint& operator=(IntegrationPoint&&) noexcept = delete;
        ~IntegrationPoint() = default;
    };

    vector<IntegrationPoint> int_pt;

public:
    Section3D(unsigned,        // tag
              unsigned,        // material tag
              double = 0.,     // area
              vec&& = {0., 0.} // eccentricity
    );

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
