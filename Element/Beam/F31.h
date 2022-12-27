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
 * @class F31
 * @brief The F31 class.
 *
 * Reference:
 *   1. https://doi.org/10.1016/0045-7949(95)00103-N
 *   2. https://doi.org/10.1061/(ASCE)0733-9445(1997)123:7(958)
 *
 * @author tlc
 * @date 12/05/2022
 * @version 1.0.0
 * @file F31.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef F31_H
#define F31_H

#include <Element/SectionElement.h>
#include <Element/Utility/B2DC.h>

class F31 final : public SectionElement3D {
    struct IntegrationPoint final {
        double coor, weight;
        unique_ptr<Section> b_section;
        mat strain_mat;
        IntegrationPoint(double, double, unique_ptr<Section>&&);
    };

    static constexpr unsigned b_node = 2, b_dof = 6, b_size = b_dof * b_node;

    static const span b_span;

    const unsigned int_pt_num, orientation_tag;

    const double length = 0., torsion_stiff = 0.;

    vector<IntegrationPoint> int_pt;

    unique_ptr<Orientation> b_trans;

    mat initial_local_flexibility;
    mat current_local_flexibility, trial_local_flexibility;
    vec current_local_deformation, trial_local_deformation;
    vec current_local_resistance, trial_local_resistance;

public:
    F31(unsigned,     // tag
        uvec&&,       // node tag
        unsigned,     // section tag
        unsigned,     // orientation tag
        unsigned = 6, // integration points
        bool = false  // nonlinear geometry switch
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int clear_status() override;
    int commit_status() override;
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
