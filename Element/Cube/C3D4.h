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
 * @class C3D4
 * @brief The C3D4 class defines C3D4 elements with non-linear support.
 *
 * C3D4 element is a four--node linear tetrahedron. A single integration point at the centre is used.
 *
 * @author tlc
 * @date 01/08/2018
 * @version 0.1.0
 * @file C3D4.h
 * @addtogroup Cube
 * @ingroup Element
 * @{
 */

#ifndef C3D4_H
#define C3D4_H

#include <Element/MaterialElement.h>

class C3D4 final : public MaterialElement3D {
    static constexpr unsigned c_node = 4, c_dof = 3, c_size = c_dof * c_node;

    const double volume = 0.;

    unique_ptr<Material> c_material;

    mat pn_pxyz, strain_mat;

public:
    C3D4(
        unsigned,    // tag
        uvec&&,      // node tag
        unsigned,    // material tag
        bool = false // nonlinear geometry switch
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;

#ifdef SUANPAN_VTK
    void Setup() override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
