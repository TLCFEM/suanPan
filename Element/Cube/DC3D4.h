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
 * @class DC3D4
 * @brief The DC3D4 class defines DC3D4 elements.
 *
 * @author tlc
 * @date 15/12/2020
 * @version 0.1.0
 * @file DC3D4.h
 * @addtogroup Cube
 * @ingroup Element
 * @{
 */

#ifndef DC3D4_H
#define DC3D4_H

#include <Element/MaterialElement.h>
#include <Element/Utility/PhaseField.h>

class DC3D4 final : public MaterialElement3D {
    static constexpr unsigned c_node = 4, c_dof = 4, c_size = c_dof * c_node;
    static const uvec u_dof, d_dof;

    const double release_rate;
    const double volume = 0.;

    mat n_mat, pn_mat, b_mat;

    unique_ptr<Material> c_material;

    PhaseField c_phase;

    double maximum_energy = 0.;

public:
    DC3D4(
        unsigned, // tag
        uvec&&,   // node tag
        unsigned, // material tag
        double,   // characteristic length
        double    // energy release rate
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
