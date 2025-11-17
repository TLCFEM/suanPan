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
#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#endif

enum class OutputType;

class vtkBase {
#ifdef SUANPAN_VTK
public:
    vtkBase() = default;
    virtual ~vtkBase() = default;

    [[nodiscard]] virtual vtkSmartPointer<vtkCell> Setup(const uvec&) const { return nullptr; }

    /**
     * Get elemental data for VTK output.
     * Produce a 6-by-n matrix, where n is the number of nodes of the element.
     * Each column represents the data of the corresponding node.
     * Extrapolation may be needed for certain element types.
     *
     * @return A 6-by-n matrix, where n is the number of nodes of the element.
     */
    virtual mat GetData(OutputType) { return {}; }

    virtual mat GetDeformation(double) { return {}; }
#endif
};

#endif

//! @}
