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

#ifndef VTKPARSER_H
#define VTKPARSER_H

#include <Recorder/OutputType.h>
#include <suanPan.h>

class DomainBase;

#ifdef SUANPAN_VTK

#include <vtkSmartPointer.h>

class vtkUnstructuredGrid;

int vtk_parser(const shared_ptr<DomainBase>&, istringstream&);

struct vtkInfo {
    OutputType type = OutputType::U;
    double scale = 1.;
    bool on_deformed = true;
    unsigned font_size = 8;
    int canvas_size[2] = {500, 500};
    bool save_file = false;
    string file_name;
    string title_name;
    bool colorbar = true;
    int material_type = -1;
    bool store_ptr = false;
    vtkSmartPointer<vtkUnstructuredGrid> grid_ptr;
};

vtkInfo vtk_process(istringstream&);

void vtk_setup(const vtkSmartPointer<vtkUnstructuredGrid>&, const vtkInfo&);

void vtk_save(const vtkSmartPointer<vtkUnstructuredGrid>&, const vtkInfo&);

void vtk_plot_node_quantity(const shared_ptr<DomainBase>&, vtkInfo);

void vtk_plot_element_quantity(const shared_ptr<DomainBase>&, vtkInfo);

int vtk_get_index(OutputType);

const char* vtk_get_name(OutputType);

#else

inline int vtk_parser(const shared_ptr<DomainBase>&, istringstream&) { return 0; }

#endif

#endif

//! @}
