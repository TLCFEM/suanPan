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
 * @fn SectionTester
 * @brief A SectionTester function.
 * @author tlc
 * @date 20/05/2023
 * @version 0.1.0
 * @file SectionTester.h
 * @addtogroup Utility
 * @{
 */

#ifndef SECTIONTESTER_H
#define SECTIONTESTER_H

#include <memory>

class DomainBase;

int test_section(const std::shared_ptr<DomainBase>&, std::istringstream&, unsigned);
int test_section_by_deformation_history(const std::shared_ptr<DomainBase>&, std::istringstream&);

#endif

//! @}
