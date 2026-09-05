/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @class NonlocalCP4
 * @brief The NonlocalCP4 class.
 *
 * @author tlc
 * @date 15/08/2026
 * @version 0.1.0
 * @file NonlocalCP4.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef NONLOCALCP4_H
#define NONLOCALCP4_H

#include <Element/MaterialElement.h>

class NonlocalCP4 final : public MaterialElement2D {
    struct IntegrationPoint {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat n_mat, b_mat;

        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&, const mat&);
    };

    static constexpr unsigned m_node{4u}, m_dof{3u}, m_size = m_dof * m_node;

    static const uvec u_dof, d_dof;

    const double thickness;

    std::vector<IntegrationPoint> int_pt;

    mat const_mat;

public:
    NonlocalCP4(
        unsigned, // tag
        uvec&&,   // node tag
        unsigned, // material tag
        double    // thickness
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    [[nodiscard]] std::vector<vec> record(OutputType) const override;

    void print() override;

#ifdef SUANPAN_VTK
    [[nodiscard]] vtkSmartPointer<vtkCell> GetCell() const override;

    mat GetData(OutputType) override;
    mat GetDeformation(double) override;
#endif
};

#endif

//! @}
