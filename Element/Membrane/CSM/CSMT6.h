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
 * @class CSMT6
 * @brief The CSMT6 class.
 *
 * @author tlc
 * @date 19/08/2021
 * @version 0.1.0
 * @file CSMT6.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CSMT6_H
#define CSMT6_H

#include <Element/MaterialElement.h>

class CSMT6 final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat b1, b2, b3;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
    };

    static constexpr unsigned m_node = 6, m_dof = 3, m_size = m_dof * m_node;

    static const uvec t_dof;
    static const uvec r_dof;

    const double thickness;

    const double area = 0.;

    vector<IntegrationPoint> int_pt;

public:
    CSMT6(
        unsigned,    // tag
        uvec&&,      // node tag
        unsigned,    // material tag
        double = 1., // thickness
        double = -1. // length
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
