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
 * @class CP6
 * @brief The CP6 class defines CPS6 CPE6 elements.
 * @author tlc
 * @date 26/10/2017
 * @version 0.1.1
 * @file CP6.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CP6_H
#define CP6_H

#include <Element/MaterialElement.h>

class CP6 final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat pn_pxy;
        sp_mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
    };

    static constexpr unsigned m_node = 6, m_dof = 2, m_size = m_dof * m_node;

    const double thickness; // thickness

    double area = 0.; // area

    vector<IntegrationPoint> int_pt;

public:
    CP6(
        unsigned,    // tag
        uvec&&,      // node tag
        unsigned,    // material tag
        double = 1., // thickness
        bool = false // nlgeom
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
