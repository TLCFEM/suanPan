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
 * @class T2D2S
 * @brief A T2D2S class
 *
 * The T2D2S class handles both linear and nonlinear problems by using a optional corotational transformation.
 *
 * @author tlc
 * @date 28/06/2018
 * @version 0.2.0
 * @file T2D2S.h
 * @addtogroup Truss
 * @ingroup Element
 * @{
 */

#ifndef T2D2S_H
#define T2D2S_H

#include <Element/SectionElement.h>

class Orientation;

class T2D2S final : public SectionElement1D {
    static constexpr unsigned t_node = 2, t_dof = 2, t_size = t_dof * t_node;

    const double length = 0.; // length of the element

    unique_ptr<Section> t_section;   // section model
    unique_ptr<Orientation> t_trans; // transformation

    const bool log_strain; // flag to indicate if to use log strain
public:
    T2D2S(
        unsigned,     // tag
        uvec&&,       // node tag
        unsigned,     // section tag
        bool = false, // nonlinear geometry switch
        bool = true   // log strain switch
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    std::vector<vec> record(OutputType) override;

    void print() override;

#ifdef SUANPAN_VTK
    void Setup() override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
