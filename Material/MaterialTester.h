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
/**
 * @fn MaterialTester
 * @brief A MaterialTester function.
 * @author tlc
 * @date 04/01/2023
 * @version 0.3.0
 * @file MaterialTester.h
 * @addtogroup Utility
 * @{
 */

#ifndef MATERIALTESTER_H
#define MATERIALTESTER_H

#include <memory>

class DomainBase;

int test_material(const std::shared_ptr<DomainBase>&, std::istringstream&, unsigned);
int test_material_with_base3d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_load(const std::shared_ptr<DomainBase>&, std::istringstream&, unsigned);
int test_material_by_load_with_base3d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_strain_history(const std::shared_ptr<DomainBase>&, std::istringstream&);
int test_material_by_stress_history(const std::shared_ptr<DomainBase>&, std::istringstream&);

#endif

//! @}
