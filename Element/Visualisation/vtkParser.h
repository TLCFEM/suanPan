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

class vtkUnstructuredGrid;

int vtk_parser(const shared_ptr<DomainBase>&, std::istringstream&);

struct vtkInfo {
    bool color_bar = true;
    bool multi_block = false;
    bool save_file = false;
    double scale = 0.;
    int canvas_size[2] = {500, 500};
    OutputType record_type = OutputType::U;
    OutputType display_type = OutputType::U;
    std::string category;
    std::string file_name;
    std::string title_name;
    unsigned font_size = 8;
};

void vtk_cell_plot(const shared_ptr<DomainBase>&, vtkInfo);

#else

inline int vtk_parser(const shared_ptr<DomainBase>&, std::istringstream&) {
    suanpan_warning("Visualisation related functionalities are not available as the current build is not compiled with the VTK support.\n");
    return 0;
}

#endif

#endif

//! @}
