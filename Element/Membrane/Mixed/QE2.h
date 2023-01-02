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
 * @class QE2
 * @brief A QE2 class.
 *
 * Reference:
 *   1. A quadrilateral mixed finite element with two enhanced strain modes [10.1002/nme.1620381102](https://doi.org/10.1002/nme.1620381102)
 *
 * @author tlc
 * @date 13/09/2018
 * @version 0.3.0
 * @file QE2.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef QE2_H
#define QE2_H

#include <Element/MaterialElement.h>

class QE2 final : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        const double factor;
        unique_ptr<Material> m_material;
        mat P, A, B, BI;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
    };

    static constexpr unsigned m_node = 4, m_dof = 2, m_size = m_dof * m_node;

    static mat mapping;

    const double thickness;

    vector<IntegrationPoint> int_pt;

    const mat mat_stiffness, iso_mapping;

    mat HT, HIL, HILI, L, LI; // constant matrices

    mat initial_qtitt, trial_qtitt, current_qtitt; // eq. 65
    vec trial_qtifi, current_qtifi;                // eq. 65
    vec trial_lambda, current_lambda;              // enhanced strain
    vec trial_alpha, current_alpha;                // strain
    vec trial_beta, current_beta;                  // stress

    vec pre_disp;

    static vec form_stress_mode(double, double);

public:
    QE2(unsigned,   // tag
        uvec&&,     // node tags
        unsigned,   // material tags
        double = 1. // thickness
        );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int clear_status() override;
    int commit_status() override;
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
