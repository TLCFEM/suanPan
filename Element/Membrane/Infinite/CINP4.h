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
 * @class CINP4
 * @brief The CINP4 class handles CINPS4 and CINPE4 elements.
 *
 * reference:
 *   1. https://doi.org/10.1016/0045-7949(84)90019-1
 *
 * @author tlc
 * @date 31/07/2018
 * @version 0.1.0
 * @file CINP4.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CINP4_H
#define CINP4_H

#include <Element/MaterialElement.h>

class CINP4 final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat pn_pxy, strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
    };

    static constexpr unsigned m_node = 4, m_dof = 2, m_size = m_dof * m_node;

    const double thickness;

    std::vector<IntegrationPoint> int_pt;

    static mat compute_mapping(const vec&);
    static mat compute_n(const vec&);
    static mat compute_dn(const vec&);
    static void stack_stiffness(mat&, const mat&, const mat&, double);

public:
    CINP4(
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
