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
 * @class B21H
 * @brief The B21H class.
 * @author tlc
 * @date 11/10/2017
 * @version 0.1.1
 * @file B21H.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef B21H_H
#define B21H_H

#include <Element/SectionElement.h>
#include <Element/Utility/B2DC.h>

class B21H final : public SectionElement2D {
    struct IntegrationPoint final {
        double coor, weight;
        unique_ptr<Section> b_section;
        mat strain_mat;
        IntegrationPoint(double, double, unique_ptr<Section>&&);
    };

    static constexpr unsigned b_node = 2, b_dof = 3, b_size = b_dof * b_node;

    const double hinge_length;

    const double length = 0.;

    vector<IntegrationPoint> int_pt, elastic_int_pt;

    unique_ptr<Orientation> b_trans;

    mat elastic_local_stiffness;

public:
    B21H(unsigned,    // tag
         uvec&&,      // node tags
         unsigned,    // section tag
         double = .2, // hinge length
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
