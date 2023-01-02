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
 * @class SGCMS
 * @brief A SGCMS shell class.
 *
 * @author tlc
 * @date 09/10/2019
 * @version 0.1.0
 * @file SGCMS.h
 * @addtogroup Shell
 * @ingroup Element
 * @{
 */

#ifndef SGCMS_H
#define SGCMS_H

#include <Element/Shell/ShellBase.h>

class SGCMS final : public ShellBase {
    struct IntegrationPoint final {
        struct SectionIntegrationPoint final {
            const double eccentricity, factor;
            unique_ptr<Material> s_material;
            SectionIntegrationPoint(double, double, unique_ptr<Material>&&);
            SectionIntegrationPoint(const SectionIntegrationPoint&);
            SectionIntegrationPoint(SectionIntegrationPoint&&) noexcept = default;
            SectionIntegrationPoint& operator=(const SectionIntegrationPoint&) = delete;
            SectionIntegrationPoint& operator=(SectionIntegrationPoint&&) noexcept = delete;
            ~SectionIntegrationPoint() = default;
        };

        vec coor;
        mat BM, BP;
        vector<SectionIntegrationPoint> sec_int_pt;
        explicit IntegrationPoint(vec&&);
    };

    static constexpr unsigned s_node = 4, s_dof = 6, s_size = s_dof * s_node;

    static const mat mapping;

    const double thickness;

    vector<IntegrationPoint> int_pt;

    static field<mat> form_plate_transformation(const mat&);
    static mat form_drilling_n(const vec&, const vec&);
    static mat form_drilling_dn(const vec&, const vec&);
    static mat form_displacement_dn(const mat&, const mat&);

public:
    SGCMS(unsigned,   // tag
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

#ifdef SUANPAN_VTK
    void Setup() override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
