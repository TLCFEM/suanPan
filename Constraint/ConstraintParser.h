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

#include <suanPan.h>

int create_new_constraint(const shared_ptr<DomainBase>&, istringstream&);

void new_embed(unique_ptr<Constraint>&, istringstream&, unsigned);
void new_fixedlength(unique_ptr<Constraint>&, istringstream&, unsigned);
void new_minimumgap(unique_ptr<Constraint>&, istringstream&, unsigned);
void new_maximumgap(unique_ptr<Constraint>&, istringstream&, unsigned);
void new_sleeve(unique_ptr<Constraint>&, istringstream&, unsigned);
void new_mpc(unique_ptr<Constraint>&, istringstream&);
void new_particlecollision2d(unique_ptr<Constraint>&, istringstream&);
void new_particlecollision3d(unique_ptr<Constraint>&, istringstream&);
void new_rigidwall(unique_ptr<Constraint>&, istringstream&, bool, bool);
void new_bc(unique_ptr<Constraint>&, istringstream&, bool, bool);
void new_nodeline(unique_ptr<Constraint>&, istringstream&);
void new_nodefacet(unique_ptr<Constraint>&, istringstream&);

#endif
