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
 * @class Allman
 * @brief The Allman class.
 * @author tlc
 * @date 14/04/2019
 * @version 0.1.0
 * @file Allman.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef ALLMAN_H
#define ALLMAN_H

#include <Element/MaterialElement.h>

class Allman final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
    };

    static constexpr unsigned m_node = 3, m_dof = 3, m_size = m_dof * m_node;

    const double thickness; // thickness

    double area = 0.; // area

    std::vector<IntegrationPoint> int_pt;

    static mat form_coor(const mat&);
    static field<mat> form_transform(const mat&);

public:
    Allman(
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
