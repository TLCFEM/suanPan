/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class NMB21
 * @brief The NMB21 class.
 * @author tlc
 * @date 26/11/2021
 * @version 0.1.0
 * @file NMB21.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef NMB21_H
#define NMB21_H

#include <Element/SectionElement.h>
#include <Element/Utility/Orientation.h>

class NMB21 : public SectionNMElement2D {
    static constexpr unsigned b_node = 2, b_dof = 3, b_size = b_dof * b_node;

protected:
    const double length = 0.;

    unique_ptr<Orientation> b_trans;
    unique_ptr<Section> b_section;

public:
    NMB21(unsigned,    // tag
          uvec&&,      // node tags
          unsigned,    // section tag
          bool = false // nonlinear geometry switch
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
