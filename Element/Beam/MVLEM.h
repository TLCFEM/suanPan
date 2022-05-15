/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class MVLEM
 * @brief The MVLEM class.
 * @author tlc
 * @date 14/09/2017
 * @version 0.1.1
 * @file MVLEM.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef MVLEM_H
#define MVLEM_H

#include <Element/MaterialElement.h>

class MVLEM final : public MaterialElement1D {
    struct Fibre final {
        double eccentricity = 0., width, height, c_area, s_area;
        unique_ptr<Material> c_material, s_material;
        Fibre(double, double, double);
    };

    static constexpr unsigned b_node = 2, b_dof = 3, b_size = b_dof * b_node;

    double shear_height;
    double length = 0.;
    double shear_height_a = 0.;
    double shear_height_b = 0.;
    double total_area = 0.;

    mat trans_mat;

    vector<Fibre> axial_spring;

    const unsigned shear_spring_tag;

    unique_ptr<Material> shear_spring;

public:
    MVLEM(unsigned,              // tag
          uvec&&,                // node tag
          const vector<double>&, // width
          const vector<double>&, // thickness
          const vector<double>&, // reinforcement ratio
          uvec&&,                // concrete material tags
          uvec&&,                // steel material tags
          unsigned,              // shear spring tag
          double                 // shear spring height
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
