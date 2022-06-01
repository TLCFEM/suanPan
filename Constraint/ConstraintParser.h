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

#ifndef CONSTRAINTPARSER_H
#define CONSTRAINTPARSER_H

#include <memory>

class DomainBase;

int create_new_constraint(const std::shared_ptr<DomainBase>&, std::istringstream&);

int create_new_bc(const std::shared_ptr<DomainBase>&, std::istringstream&, bool);
int create_new_groupbc(const std::shared_ptr<DomainBase>&, std::istringstream&, bool);
int create_new_fixedlength(const std::shared_ptr<DomainBase>&, std::istringstream&, unsigned);
int create_new_mpc(const std::shared_ptr<DomainBase>&, std::istringstream&);
int create_new_particlecollision2d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int create_new_particlecollision3d(const std::shared_ptr<DomainBase>&, std::istringstream&);
int create_new_rigidwall(const std::shared_ptr<DomainBase>&, std::istringstream&, bool, bool);

int create_new_criterion(const std::shared_ptr<DomainBase>&, std::istringstream&);

#endif
