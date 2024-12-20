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
 * @class Section2D
 * @brief A Section2D class.
 * @author tlc
 * @date 13/03/2019
 * @version 0.2.0
 * @file Section2D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef SECTION2D_H
#define SECTION2D_H

#include <Section/Section.h>
#include <Material/Material.h>
#include <Toolbox/ResourceHolder.h>

using std::vector;

class Section2D : public Section {
protected:
    struct IntegrationPoint {
        double coor, weight;
        ResourceHolder<Material> s_material;
        IntegrationPoint(double, double, unique_ptr<Material>&&);
    };

    vector<IntegrationPoint> int_pt;

    void initialize_stiffness();

public:
    Section2D(
        unsigned,    // tag
        unsigned,    // material tag
        double = 0., // area
        double = 0.  // eccentricity
    );

    void set_characteristic_length(double) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
