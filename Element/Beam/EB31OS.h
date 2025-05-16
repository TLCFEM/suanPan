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
 * @class EB31OS
 * @brief The EB31OS class.
 *
 * Elastic 3D beam element with additional warping DOFs.
 *
 * @author tlc
 * @date 10/09/2023
 * @version 0.1.0
 * @file EB31OS.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef EB31OS_H
#define EB31OS_H

#include <Element/SectionElement.h>
#include <Element/Utility/Orientation.h>

class EB31OS final : public SectionOSElement3D {
    static constexpr unsigned b_node = 2u, b_dof = 7u, b_size = b_dof * b_node;

    const unsigned orientation_tag;

    const double length = 0.;

    unique_ptr<Orientation> b_trans;

    const vec property; // [E, G, A, IZ, IY, J, IW]

    mat local_stiff;

public:
    EB31OS(
        unsigned,    // tag
        uvec&&,      // node tags
        vec&&,       // properties
        unsigned,    // orientation tag
        bool = false // nonlinear geometry switch
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    std::vector<vec> record(OutputType) override;

    void print() override;

#ifdef SUANPAN_VTK
    void Setup() override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
