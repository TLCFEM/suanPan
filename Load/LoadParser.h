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

#ifndef LOADPARSER_H
#define LOADPARSER_H

#include <memory>

class DomainBase;

int create_new_acceleration(const std::shared_ptr<DomainBase>&, std::istringstream&);
int create_new_amplitude(const std::shared_ptr<DomainBase>&, std::istringstream&);
int create_new_bodyforce(const std::shared_ptr<DomainBase>&, std::istringstream&, bool);
int create_new_cload(const std::shared_ptr<DomainBase>&, std::istringstream&, bool = false);
int create_new_lineudl(const std::shared_ptr<DomainBase>&, std::istringstream&, unsigned);
int create_new_displacement(const std::shared_ptr<DomainBase>&, std::istringstream&, bool = false);
int create_new_supportmotion(const std::shared_ptr<DomainBase>&, std::istringstream&, unsigned);

#endif
