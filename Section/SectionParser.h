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

#ifndef SECTIONPARSER_H
#define SECTIONPARSER_H

#include <suanPan.h>

using std::istringstream;

class Section;
class DomainBase;

int create_new_section(const shared_ptr<DomainBase>&, istringstream&);

// 1D
void new_circle1d(unique_ptr<Section>&, istringstream&);
void new_fibre1d(unique_ptr<Section>&, istringstream&);
void new_rectangle1d(unique_ptr<Section>&, istringstream&);
void new_trusssection(unique_ptr<Section>&, istringstream&);

// 2D
void new_bar2d(unique_ptr<Section>&, istringstream&);
void new_box2d(unique_ptr<Section>&, istringstream&);
void new_circle2d(unique_ptr<Section>&, istringstream&);
void new_circularhollow2D(unique_ptr<Section>&, istringstream&);
void new_fibre2d(unique_ptr<Section>&, istringstream&);
void new_hsection2d(unique_ptr<Section>&, istringstream&);
void new_isection2d(unique_ptr<Section>&, istringstream&);
void new_rectangle2d(unique_ptr<Section>&, istringstream&);
void new_tsection2d(unique_ptr<Section>&, istringstream&);

void new_nz2d(unique_ptr<Section>&, istringstream&);
void new_us2d(unique_ptr<Section>&, istringstream&);
void new_eu2d(unique_ptr<Section>&, istringstream&);

// 3D
void new_bar3d(unique_ptr<Section>&, istringstream&);
void new_box3d(unique_ptr<Section>&, istringstream&);
void new_circle3d(unique_ptr<Section>&, istringstream&);
void new_circularhollow3D(unique_ptr<Section>&, istringstream&);
void new_fibre3d(unique_ptr<Section>&, istringstream&);
void new_isection3d(unique_ptr<Section>&, istringstream&);
void new_rectangle3d(unique_ptr<Section>&, istringstream&);
void new_tsection3d(unique_ptr<Section>&, istringstream&);

void new_nz3d(unique_ptr<Section>&, istringstream&);
void new_us3d(unique_ptr<Section>&, istringstream&);
void new_eu3d(unique_ptr<Section>&, istringstream&);

// NM
void new_nm2d1(unique_ptr<Section>&, istringstream&);
void new_nm2d2(unique_ptr<Section>&, istringstream&);
void new_nm3d1(unique_ptr<Section>&, istringstream&);
void new_nm3d2(unique_ptr<Section>&, istringstream&);

vec euisection(const string&);
vec nzchsection(const string&);
vec nzisection(const string&);
vec nzrhsection(const string&);
vec nzshsection(const string&);
vec usisection(const string&);
vec ustsection(const string&);

#endif
