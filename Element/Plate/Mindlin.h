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
 * @class Mindlin
 * @brief A Mindlin plate class.
 *
 * @author tlc
 * @date 04/07/2018
 * @version 0.1.0
 * @file Mindlin.h
 * @addtogroup Plate
 * @ingroup Element
 * @{
 */

#ifndef MINDLIN_H
#define MINDLIN_H

#include <Element/MaterialElement.h>

class Mindlin final : public MaterialElement2D {
    struct IntegrationPoint final {
        struct SectionIntegrationPoint final {
            const double eccentricity, factor;
            unique_ptr<Material> p_material;
            SectionIntegrationPoint(double, double, unique_ptr<Material>&&);
            SectionIntegrationPoint(const SectionIntegrationPoint&);
        };

        vec coor;
        mat strain_mat;
        vector<SectionIntegrationPoint> sec_int_pt;
        explicit IntegrationPoint(vec&&);
    };

    static constexpr unsigned p_node = 4, p_dof = 3, p_size = p_dof * p_node;

    const double thickness;
    const unsigned num_section_ip;

    vector<IntegrationPoint> int_pt;

    mat penalty_stiffness;

public:
    Mindlin(
        unsigned,    // element tag
        uvec&&,      // node tag
        unsigned,    // material tag
        double,      // thickness
        unsigned = 2 // integration points along thickness
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
