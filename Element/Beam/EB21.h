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
 * @class EB21
 * @brief The ElasticB21 class.
 * @author tlc
 * @date 14/09/2017
 * @version 0.1.1
 * @file EB21.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef ELASTICB21_H
#define ELASTICB21_H

#include <Element/MaterialElement.h>
#include <Element/Utility/Orientation.h>

class EB21 final : public MaterialElement1D {
    static constexpr unsigned b_node = 2, b_dof = 3, b_size = b_dof * b_node;

    const double area, moment_inertia;

    unique_ptr<Material> b_material;
    unique_ptr<Orientation> b_trans;

    mat local_stiff;

public:
    EB21(
        unsigned,    // tag
        uvec&&,      // node tag
        double,      // area
        double,      // moment of inertia
        unsigned,    // material tags
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
