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
 * @class DCP3
 * @brief The DCP3 class.
 * @author tlc
 * @date 08/12/2020
 * @version 0.1.0
 * @file DCP3.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef DCP3_H
#define DCP3_H

#include <Element/MaterialElement.h>
#include <Element/Utility/PhaseField.h>

class DCP3 final : public MaterialElement2D {
    static constexpr unsigned m_node = 3, m_dof = 3, m_size = m_dof * m_node;
    static const uvec u_dof, d_dof;

    const double release_rate;
    const double thickness; // thickness
    const double area = 0.; // area

    mat n_mat, pn_mat, b_mat;

    unique_ptr<Material> m_material; // store material model

    PhaseField m_phase;

    double maximum_energy = 0.;

public:
    DCP3(unsigned,   // tag
         uvec&&,     // node tag
         unsigned,   // material tag
         double,     // characteristic length
         double,     // energy release rate
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
    mat GetData(OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
