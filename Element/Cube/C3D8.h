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
 * @class C3D8
 * @brief The C3D8 class defines C3D8 C3D8R elements.
 * @author tlc
 * @date 13/07/2018
 * @version 0.2.0
 * @file C3D8.h
 * @addtogroup Cube
 * @ingroup Element
 * @{
 */

#ifndef C3D8_H
#define C3D8_H

#include <Element/MaterialElement.h>

class C3D8 final : public MaterialElement3D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> c_material;
        mat pn_pxyz;
        sp_mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
    };

    static constexpr unsigned c_node = 8, c_dof = 3, c_size = c_dof * c_node;

    static const field<vec> h_mode;

    const char int_scheme;

    const bool hourglass_control;

    mat hourglass;

    std::vector<IntegrationPoint> int_pt;

public:
    C3D8(
        unsigned,    // tag
        uvec&&,      // node tag
        unsigned,    // material tag
        char = 'I',  // reduced integration
        bool = false // nonlinear geometry switch
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    [[nodiscard]] mat compute_shape_function(const mat&, unsigned) const override;

    std::vector<vec> record(OutputType) override;

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
