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
 * @class C3D8I
 * @brief The C3D8I class defines C3D8I C3D8IR elements.
 * @author tlc
 * @date 26/07/2018
 * @version 0.1.0
 * @file C3D8I.h
 * @addtogroup Cube
 * @ingroup Element
 * @{
 */

#ifndef C3D8I_H
#define C3D8I_H

#include <Element/MaterialElement.h>

class C3D8I final : public MaterialElement3D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> c_material;
        mat pn_pxyz, B1, B2;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
    };

    static constexpr unsigned c_node = 8, c_dof = 3, c_size = c_dof * c_node;

    vector<IntegrationPoint> int_pt;

public:
    C3D8I(unsigned, // tag
          uvec&&,   // node tag
          unsigned  // material tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    [[nodiscard]] mat compute_shape_function(const mat&, unsigned) const override;

    vector<vec> record(OutputType) override;

    void print() override;

#ifdef SUANPAN_VTK
    void Setup() override;
    mat GetData(OutputType) override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
