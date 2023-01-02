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
 * @class PatchQuad
 * @brief The PatchQuad class handles CPS4, CPE4, CPS4R and CPE4R elements. It is a
 * four node constant strain membrane element with optional reduced integration
 * for both plane stress and plane strain problems and optional switch for TL
 * nonlinear geometry formulation.
 *
 * @author tlc
 * @date 15/11/2020
 * @version 0.1.0
 * @file PatchQuad.h
 * @addtogroup Patch
 * @ingroup Element
 * @{
 */

#ifndef PATCHQUAD_H
#define PATCHQUAD_H

#include <Element/Patch/Patch.h>

class PatchQuad final : public MaterialPatch2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        sp_mat strain_mat;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
    };

    static constexpr unsigned m_dof = 2;

    const unsigned m_node, m_size;

    const double thickness;

    vector<IntegrationPoint> int_pt;

public:
    PatchQuad(unsigned,   // tag
              vec&&,      // knot x
              vec&&,      // knot y
              uvec&&,     // node tag
              unsigned,   // material tag
              double = 1. // thickness
        );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
