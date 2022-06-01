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

#ifndef MATERIALPARSER_H
#define MATERIALPARSER_H

#include <memory>

class DomainBase;

int create_new_material(const std::shared_ptr<DomainBase>&, std::istringstream&);

int test_material1d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material2d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material3d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_with_base3d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_load1d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_load2d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_load3d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_load_with_base3d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_strain_history(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_stress_history(const std::shared_ptr<DomainBase>&, std::istringstream&);

#endif
