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
 * @class DC3D8
 * @brief The DC3D8 class.
 * @author tlc
 * @date 16/12/2020
 * @version 0.1.0
 * @file DC3D8.h
 * @addtogroup Cube
 * @ingroup Element
 * @{
 */

#ifndef DC3D8_H
#define DC3D8_H

#include <Element/MaterialElement.h>
#include <Element/Utility/PhaseField.h>

class DC3D8 final : public MaterialElement3D {
    struct IntegrationPoint final : PhaseField {
        vec coor;
        double weight;
        double maximum_energy = 0.;
        unique_ptr<Material> c_material;
        mat n_mat, pn_mat;
        sp_mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&, mat&&);
    };

    static constexpr unsigned c_node = 8, c_dof = 4, c_size = c_dof * c_node;

    static const uvec u_dof, d_dof;

    const double release_rate;

    vector<IntegrationPoint> int_pt;

public:
    DC3D8(unsigned, // tag
          uvec&&,   // node tag
          unsigned, // material tag
          double,   // characteristic length
          double    // release rate
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
    mat GetData(OutputType) override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
