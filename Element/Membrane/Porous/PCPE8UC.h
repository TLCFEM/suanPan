/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class PCPE8UC
 * @brief The PCPE8 class handles CPE8 elements with pore pressure.
 *
 * @author tlc
 * @date 23/11/2021
 * @version 0.1.0
 * @file PCPE8UC.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef PCPE8UC_H
#define PCPE8UC_H

#include <Element/MaterialElement.h>

class PCPE8UC final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
    };

    static constexpr unsigned m_node = 8, m_dof = 2, m_size = m_dof * m_node;

    const double alpha;
    const double porosity;

    double q = 0.;

    mat meta_k;

    vector<IntegrationPoint> int_pt;

public:
    PCPE8UC(
        unsigned, // tag
        uvec&&,   // node tag
        unsigned, // solid material tag
        unsigned, // fluid material tag
        double,   // alpha
        double    // porosity
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
