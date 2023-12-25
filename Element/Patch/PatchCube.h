/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class PatchCube
 * @brief The PatchCube class.
 * @author tlc
 * @date 13/11/2020
 * @version 0.2.0
 * @file PatchCube.h
 * @addtogroup Patch
 * @ingroup Element
 * @{
 */

#ifndef PATCHCUBE_H
#define PATCHCUBE_H

#include <Element/Patch/Patch.h>

class PatchCube final : public MaterialPatch3D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> c_material;
        sp_mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
    };

    static constexpr unsigned c_dof = 3;

    const unsigned c_node;

    const unsigned c_size = c_dof * c_node;

    vector<IntegrationPoint> int_pt;

public:
    PatchCube(
        unsigned, // tag
        vec&&,    // knot x
        vec&&,    // knot y
        vec&&,    // knot z
        uvec&&,   // node tag
        unsigned  // material tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
