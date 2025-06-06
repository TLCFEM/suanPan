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

#ifndef VTKPARSER_H
#define VTKPARSER_H

#include <Recorder/OutputType.h>
#include <suanPan.h>

class DomainBase;

#ifdef SUANPAN_VTK

#include <vtkSmartPointer.h>

class vtkUnstructuredGrid;

int vtk_parser(const shared_ptr<DomainBase>&, std::istringstream&);

struct vtkInfo {
    OutputType type = OutputType::U;
    double scale = 1.;
    bool on_deformed = true;
    unsigned font_size = 8;
    int canvas_size[2] = {500, 500};
    bool save_file = false;
    std::string file_name;
    std::string title_name;
    bool colorbar = true;
    int material_type = -1;
    bool store_ptr = false;
    vtkSmartPointer<vtkUnstructuredGrid> grid_ptr;
};

vtkInfo vtk_process(std::istringstream&);

void vtk_setup(const vtkSmartPointer<vtkUnstructuredGrid>&, const vtkInfo&);

void vtk_save(vtkSmartPointer<vtkUnstructuredGrid>&&, vtkInfo);

void vtk_plot_node_quantity(const shared_ptr<DomainBase>&, vtkInfo);

void vtk_plot_element_quantity(const shared_ptr<DomainBase>&, vtkInfo);

#else

inline int vtk_parser(const shared_ptr<DomainBase>&, std::istringstream&) {
    suanpan_warning("Visualisation related functionalities are not available as the current build is not compiled with the VTK support.\n");
    return 0;
}

#endif

#endif

//! @}
