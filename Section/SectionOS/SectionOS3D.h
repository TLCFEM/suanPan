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
 * @class SectionOS3D
 * @brief A SectionOS3D class.
 *
 * Reference:
 *   1. Distributed plasticity analysis of steel building structural systems
 *
 * @author tlc
 * @date 15/09/2023
 * @version 0.1.0
 * @file SectionOS3D.h
 * @addtogroup Section-OS
 * @ingroup Section
 * @{
 */

#ifndef SECTIONOS3D_H
#define SECTIONOS3D_H

#include <Section/Section.h>
#include <Material/Material.h>
#include <Toolbox/ResourceHolder.h>

using std::vector;

class SectionOS3D : public Section {
    static const mat weighing_mat;

protected:
    struct IntegrationPoint {
        const double omega, py, pz, weight;
        ResourceHolder<Material> s_material;
        IntegrationPoint(double, double, double, double, unique_ptr<Material>&&);
    };

    vector<IntegrationPoint> int_pt;

public:
    SectionOS3D(unsigned,        // tag
                unsigned,        // material tag
                double = 0.,     // area
                vec&& = {0., 0.} // eccentricity
    );

    void set_characteristic_length(double) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;
};

#endif

//! @}
