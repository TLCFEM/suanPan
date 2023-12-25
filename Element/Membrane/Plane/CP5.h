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
 * @class CP5
 * @brief The CP5 class handles CPS5, CPE5, CPS5R and CPE5R elements. It is a
 * four node constant strain membrane element with optional reduced integration
 * for both plane stress and plane strain problems.
 *
 * @author tlc
 * @date 08/02/2020
 * @version 0.1.0
 * @file CP5.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CP5_H
#define CP5_H

#include <Element/MaterialElement.h>

class CP5 final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat pn_pxy;
        sp_mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
    };

    static constexpr unsigned m_node = 5, m_dof = 2, m_size = m_dof * m_node;

    const double thickness;

    vector<IntegrationPoint> int_pt;

public:
    CP5(
        unsigned,    // tag
        uvec&&,      // node tag
        unsigned,    // material tag
        double = 1., // thickness
        bool = false // nonlinear geometry switch
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    [[nodiscard]] mat compute_shape_function(const mat&, unsigned) const override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
