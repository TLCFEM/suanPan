/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @class NonlocalC3D8
 * @brief The NonlocalC3D8 class.
 * @author tlc
 * @date 16/12/2020
 * @version 0.1.0
 * @file NonlocalC3D8.h
 * @addtogroup Cube
 * @ingroup Element
 * @{
 */

#ifndef NONLOCALC3D8_H
#define NONLOCALC3D8_H

#include <Element/MaterialElement.h>

class NonlocalC3D8 final : public MaterialElement3D {
    struct IntegrationPoint {
        vec coor;
        double weight;
        unique_ptr<Material> c_material;
        sp_mat strain_mat;

        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, const mat&, const mat&);
    };

    static constexpr unsigned c_node{8u}, c_dof{4u}, c_size = c_dof * c_node;

    static const uvec u_dof, d_dof;

    std::vector<IntegrationPoint> int_pt;

    mat const_mat;

public:
    NonlocalC3D8(
        unsigned, // tag
        uvec&&,   // node tag
        unsigned  // material tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    [[nodiscard]] std::vector<vec> record(OutputType) const override;

    void print() override;

#ifdef SUANPAN_VTK
    [[nodiscard]] vtkSmartPointer<vtkCell> GetCell() const override;

    mat GetData(OutputType) override;
    mat GetDeformation(double) override;
#endif
};

#endif

//! @}
