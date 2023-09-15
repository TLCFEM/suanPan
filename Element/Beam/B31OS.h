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
 * @class B31OS
 * @brief The B31OS class.
 * @author tlc
 * @date 15/09/2023
 * @version 0.1.0
 * @file B31OS.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef B31OS_H
#define B31OS_H

#include <Element/SectionElement.h>
#include <Element/Utility/Orientation.h>

class B31OS final : public SectionOSElement3D {
    struct IntegrationPoint final {
        double coor, weight;
        unique_ptr<Section> b_section;
        sp_mat strain_mat;
        IntegrationPoint(double, double, double, unique_ptr<Section>&&);
    };

    static constexpr unsigned b_node = 2, b_dof = 7, b_size = b_dof * b_node;

    const unsigned orientation_tag, int_pt_num;

    const double length{0.};

    vector<IntegrationPoint> int_pt;

    unique_ptr<Orientation> b_trans;

public:
    B31OS(unsigned,     // tag
          uvec&&,       // node tags
          unsigned,     // section tag
          unsigned,     // orientation tag
          unsigned = 6, // integration points
          bool = false  // nonlinear geometry switch
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
