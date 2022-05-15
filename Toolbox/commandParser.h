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

#ifndef COMMANDPARSER_H
#define COMMANDPARSER_H

#include <suanPan.h>

class Bead;
class DomainBase;

int process_command(const shared_ptr<Bead>&, istringstream&);

int process_file(const shared_ptr<Bead>&, const char*);
int process_file(const shared_ptr<Bead>&, istringstream&);

int run_example();

int create_new_domain(const shared_ptr<Bead>&, istringstream&);

int disable_object(const shared_ptr<Bead>&, istringstream&);
int enable_object(const shared_ptr<Bead>&, istringstream&);
int erase_object(const shared_ptr<Bead>&, istringstream&);

int save_object(const shared_ptr<DomainBase>&, istringstream&);
int list_object(const shared_ptr<DomainBase>&, istringstream&);
int suspend_object(const shared_ptr<DomainBase>&, istringstream&);
int protect_object(const shared_ptr<DomainBase>&, istringstream&);

int create_new_acceleration(const shared_ptr<DomainBase>&, istringstream&);
int create_new_amplitude(const shared_ptr<DomainBase>&, istringstream&);
int create_new_bc(const shared_ptr<DomainBase>&, istringstream&, bool);
int create_new_bodyforce(const shared_ptr<DomainBase>&, istringstream&, bool);
int create_new_groupbc(const shared_ptr<DomainBase>&, istringstream&, bool);
int create_new_cload(const shared_ptr<DomainBase>&, istringstream&, bool = false);
int create_new_lineudl(const shared_ptr<DomainBase>&, istringstream&, unsigned);
int create_new_converger(const shared_ptr<DomainBase>&, istringstream&);
int create_new_criterion(const shared_ptr<DomainBase>&, istringstream&);
int create_new_displacement(const shared_ptr<DomainBase>&, istringstream&, bool = false);
int create_new_external_module(const shared_ptr<DomainBase>&, istringstream&);
int create_new_fixedlength(const shared_ptr<DomainBase>&, istringstream&, unsigned);
int create_new_nodegroup(const shared_ptr<DomainBase>&, istringstream&);
int create_new_elementgroup(const shared_ptr<DomainBase>&, istringstream&);
int create_new_generate(const shared_ptr<DomainBase>&, istringstream&);
int create_new_generatebyrule(const shared_ptr<DomainBase>&, istringstream&);
int create_new_generatebyplane(const shared_ptr<DomainBase>&, istringstream&);
int create_new_generatebypoint(const shared_ptr<DomainBase>&, istringstream&);
int create_new_groupgroup(const shared_ptr<DomainBase>&, istringstream&);
int create_new_initial(const shared_ptr<DomainBase>&, istringstream&);
int create_new_integrator(const shared_ptr<DomainBase>&, istringstream&);
int create_new_mass(const shared_ptr<DomainBase>&, istringstream&);
int create_new_modifier(const shared_ptr<DomainBase>&, istringstream&);
int create_new_mpc(const shared_ptr<DomainBase>&, istringstream&);
int create_new_node(const shared_ptr<DomainBase>&, istringstream&);
int create_new_orientation(const shared_ptr<DomainBase>&, istringstream&);
int create_new_recorder(const shared_ptr<DomainBase>&, istringstream&);
int create_new_rigidwall(const shared_ptr<DomainBase>&, istringstream&, bool, bool);
int create_new_particlecollision2d(const shared_ptr<DomainBase>&, istringstream&);
int create_new_particlecollision3d(const shared_ptr<DomainBase>&, istringstream&);
int create_new_plainrecorder(const shared_ptr<DomainBase>&, istringstream&);
int create_new_hdf5recorder(const shared_ptr<DomainBase>&, istringstream&);
int create_new_solver(const shared_ptr<DomainBase>&, istringstream&);
int create_new_step(const shared_ptr<DomainBase>&, istringstream&);
int create_new_supportmotion(const shared_ptr<DomainBase>&, istringstream&, unsigned);

int test_material1d(const shared_ptr<DomainBase>&, istringstream&);
int test_material2d(const shared_ptr<DomainBase>&, istringstream&);
int test_material3d(const shared_ptr<DomainBase>&, istringstream&);
int test_material_with_base3d(const shared_ptr<DomainBase>&, istringstream&);
int test_material_by_load1d(const shared_ptr<DomainBase>&, istringstream&);
int test_material_by_load2d(const shared_ptr<DomainBase>&, istringstream&);
int test_material_by_load3d(const shared_ptr<DomainBase>&, istringstream&);
int test_material_by_load_with_base3d(const shared_ptr<DomainBase>&, istringstream&);
int test_material_by_strain_history(const shared_ptr<DomainBase>&, istringstream&);
int test_material_by_stress_history(const shared_ptr<DomainBase>&, istringstream&);

int set_property(const shared_ptr<DomainBase>&, istringstream&);

int print_info(const shared_ptr<DomainBase>&, istringstream&);

int print_command();
int execute_command(istringstream&);

#endif
