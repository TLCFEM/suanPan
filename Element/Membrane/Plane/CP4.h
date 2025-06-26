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
 * @class CP4
 * @brief The CP4 class handles CPS4, CPE4, CPS4R and CPE4R elements. It is a
 * four node constant strain membrane element with optional reduced integration
 * for both plane stress and plane strain problems and optional switch for TL
 * nonlinear geometry formulation.
 *
 * References:
 *     1. [10.1016/0045-7825(84)90067-7](https://doi.org/10.1016/0045-7825(84)90067-7)
 *
 * @author tlc
 * @date 10/06/2018
 * @version 0.1.3
 * @file CP4.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CP4_H
#define CP4_H

#include <Element/MaterialElement.h>

class CP4 final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat pn_pxy;
        sp_mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
    };

    static constexpr unsigned m_node = 4, m_dof = 2, m_size = m_dof * m_node;

    static const vec h_mode;

    const double thickness, penalty;

    const bool reduced_scheme;

    std::vector<IntegrationPoint> int_pt;

    mat hourglass;

    static void stack_stiffness(mat&, const mat&, const sp_mat&, double);

public:
    CP4(
        unsigned, // tag
        uvec&&,   // node tag
        unsigned, // material tag
        double,   // thickness
        double,   // multiplier for hourglassing control
        bool,     // reduced integration
        bool      // nonlinear geometry switch
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    [[nodiscard]] mat compute_shape_function(const mat&, unsigned) const override;

    std::vector<vec> record(OutputType) override;

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
