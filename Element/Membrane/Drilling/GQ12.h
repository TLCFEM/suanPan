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
 * @class GQ12
 * @brief The GQ12 class implements the displacement based four node
 * quadrilateral drilling element proposed by Long and Xu (1994).
 *
 * Reference:
 * 1. Generalized conforming Quadrilateral Membrane Element With Vertex Rigid
 * Rotational Freedom. https://doi.org/10.1016/0045-7949(94)90356-5
 *
 * The element assumes the displacement field is compatible/conforming on
 * element boundaries in an averaged/weak sense. The element exhibits a good
 * performance.
 *
 * @author tlc
 * @date 26/01/2018
 * @version 0.1.2
 * @file GQ12.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef GQ12_H
#define GQ12_H

#include <Element/MaterialElement.h>

class GQ12 final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
    };

    static constexpr unsigned m_node = 4, m_dof = 3, m_size = m_dof * m_node;

    const double thickness;

    vector<IntegrationPoint> int_pt;

public:
    GQ12(
        unsigned,   // tag
        uvec&&,     // node tag
        unsigned,   // material tag
        double = 1. // thickness
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    [[nodiscard]] mat compute_shape_function(const mat&, unsigned) const override;

    vector<vec> record(OutputType) override;

    void print() override;

#ifdef SUANPAN_VTK
    void Setup() override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    mat GetData(OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
