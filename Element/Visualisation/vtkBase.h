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
 * @class vtkBase
 * @brief The vtkBase class.
 *
 * @author tlc
 * @date 28/07/2018
 * @version 0.1.3
 * @file vtkBase.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef VTKBASE_H
#define VTKBASE_H

#ifdef SUANPAN_VTK
#include <Recorder/OutputType.h>
#include <armadillo/armadillo>
#include <vtkCell.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

using arma::mat;
#endif

class vtkBase {
#ifdef SUANPAN_VTK

protected:
    vtkSmartPointer<vtkCell> vtk_cell;

public:
    vtkBase() = default;
    virtual ~vtkBase() = default;

    virtual void Setup();

    virtual void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType);
    virtual mat GetData(OutputType);

    virtual void SetDeformation(vtkSmartPointer<vtkPoints>&, double);

    [[nodiscard]] virtual const vtkSmartPointer<vtkCell>& GetCell() const;
#endif
};

#endif

//! @}
