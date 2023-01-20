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

#include "Orientation.h"
#include <Element/Element.h>

void Orientation::check_element_ptr() const { suanpan_assert([&] { if(element_ptr == nullptr) throw logic_error("need to set element pointer first"); }); }

Orientation::Orientation(const unsigned T, vec&& O)
    : Tag(T)
    , z_axis(std::forward<vec>(O)) {}

void Orientation::update_axis(const vec& O) { if(O.n_elem == 3) z_axis = O; }

void Orientation::set_element_ptr(const Element* E) {
    element_ptr = E;
    update_status();
}

bool Orientation::is_nlgeom() const { return false; }

double Orientation::get_length() const { return length; }

double Orientation::get_inclination() const { return inclination; }

const mat& Orientation::get_transformation() const { return direction_cosine; }

void Orientation::update_status() {
    check_element_ptr();
    update_transformation();
}

void Orientation::commit_status() {}

void Orientation::reset_status() {}

void Orientation::clear_status() {}

vec Orientation::to_local_vec(const double g) const { return to_local_vec(vec{g}); }

vec Orientation::to_global_vec(const double l) const { return to_global_vec(vec{l}); }

mat Orientation::to_global_mass_mat(const double l) const { return to_global_mass_mat(mat{l}); }

mat Orientation::to_global_geometry_mat(const double l) const { return to_global_geometry_mat(mat{l}); }

mat Orientation::to_global_stiffness_mat(const double l) const { return to_global_stiffness_mat(mat{l}); }

vec Orientation::to_local_vec(vec&& in) const { return to_local_vec(in); }

vec Orientation::to_global_vec(vec&& in) const { return to_global_vec(in); }

mat Orientation::to_global_mass_mat(mat&& in) const { return to_global_mass_mat(in); }

mat Orientation::to_global_geometry_mat(mat&& in) const { return to_global_geometry_mat(in); }

mat Orientation::to_global_stiffness_mat(mat&& in) const { return to_global_stiffness_mat(in); }

mat Orientation::to_global_mass_mat(const mat&) const { throw logic_error("not implemented.\n"); }

mat Orientation::to_global_geometry_mat(const mat&) const { throw logic_error("not applicable to linear formulation.\n"); }
