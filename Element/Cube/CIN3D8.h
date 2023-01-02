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
 * @class CIN3D8
 * @brief The CIN3D8 class defines CIN3D8 CIN3D8R elements.
 * @author tlc
 * @date 13/07/2018
 * @version 0.2.0
 * @file CIN3D8.h
 * @addtogroup Cube
 * @ingroup Element
 * @{
 */

#ifndef CIN3D8_H
#define CIN3D8_H

#include <Element/MaterialElement.h>

class CIN3D8 final : public MaterialElement3D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> c_material;
        mat pn_pxyz, strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
    };

    static constexpr unsigned c_node = 8, c_dof = 3, c_size = c_dof * c_node;

    vector<IntegrationPoint> int_pt;

    static mat compute_mapping(const vec&);
    static mat compute_n(const vec&);
    static mat compute_dn(const vec&);

public:
    CIN3D8(unsigned, // tag
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

#ifdef SUANPAN_VTK
    void Setup() override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
