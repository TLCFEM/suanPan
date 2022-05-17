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

#include <suanPan.h>

class Material;
class DomainBase;

int create_new_material(const shared_ptr<DomainBase>&, istringstream&);

void new_afc01(unique_ptr<Material>&, istringstream&);
void new_afc02(unique_ptr<Material>&, istringstream&);
void new_afc03(unique_ptr<Material>&, istringstream&);
void new_armstrongfrederick(unique_ptr<Material>&, istringstream&);
void new_armstrongfrederick1d(unique_ptr<Material>&, istringstream&);
void new_axisymmetricelastic(unique_ptr<Material>&, istringstream&);
void new_bilinear1d(unique_ptr<Material>&, istringstream&);
void new_bilinear2d(unique_ptr<Material>&, istringstream&);
void new_bilinearcc(unique_ptr<Material>&, istringstream&);
void new_bilineardp(unique_ptr<Material>&, istringstream&);
void new_bilinearelastic1d(unique_ptr<Material>&, istringstream&);
void new_bilinearhoffman(unique_ptr<Material>&, istringstream&);
void new_bilinearj2(unique_ptr<Material>&, istringstream&);
void new_bilinearmises1d(unique_ptr<Material>&, istringstream&);
void new_bilinearoo(unique_ptr<Material>&, istringstream&);
void new_bilinearpo(unique_ptr<Material>&, istringstream&);
void new_bilinearperic(unique_ptr<Material>&, istringstream&);
void new_blatzko(unique_ptr<Material>&, istringstream&);
void new_boucwen(unique_ptr<Material>&, istringstream&);
void new_bwbn(unique_ptr<Material>&, istringstream&);
void new_cdp(unique_ptr<Material>&, istringstream&);
void new_cdpm2(unique_ptr<Material>&, istringstream&, unsigned);
void new_concrete21(unique_ptr<Material>&, istringstream&);
void new_concrete22(unique_ptr<Material>&, istringstream&);
void new_concretecm(unique_ptr<Material>&, istringstream&);
void new_concreteexp(unique_ptr<Material>&, istringstream&);
void new_concretetable(unique_ptr<Material>&, istringstream&);
void new_concretetsai(unique_ptr<Material>&, istringstream&);
void new_coulombfriction(unique_ptr<Material>&, istringstream&);
void new_dafaliasmanzari(unique_ptr<Material>&, istringstream&);
void new_elastic1d(unique_ptr<Material>&, istringstream&);
void new_elastic2d(unique_ptr<Material>&, istringstream&);
void new_expcc(unique_ptr<Material>&, istringstream&);
void new_expdp(unique_ptr<Material>&, istringstream&);
void new_expgurson(unique_ptr<Material>&, istringstream&);
void new_expgurson1d(unique_ptr<Material>&, istringstream&);
void new_exphoffman(unique_ptr<Material>&, istringstream&);
void new_expj2(unique_ptr<Material>&, istringstream&);
void new_expmises1d(unique_ptr<Material>&, istringstream&);
void new_flag01(unique_ptr<Material>&, istringstream&);
void new_flag02(unique_ptr<Material>&, istringstream&);
void new_fluid(unique_ptr<Material>&, istringstream&);
void new_gap01(unique_ptr<Material>&, istringstream&);
void new_isotropicelastic3d(unique_ptr<Material>&, istringstream&);
void new_kelvin(unique_ptr<Material>&, istringstream&);
void new_lineardamage(unique_ptr<Material>&, istringstream&);
void new_maxwell(unique_ptr<Material>&, istringstream&);
void new_mooneyrivlin(unique_ptr<Material>&, istringstream&);
void new_mpf(unique_ptr<Material>&, istringstream&);
void new_multilinearoo(unique_ptr<Material>&, istringstream&);
void new_multilinearpo(unique_ptr<Material>&, istringstream&);
void new_multilinearelastic1d(unique_ptr<Material>&, istringstream&);
void new_multilinearj2(unique_ptr<Material>&, istringstream&);
void new_multilinearmises1d(unique_ptr<Material>&, istringstream&);
void new_nle1d01(unique_ptr<Material>&, istringstream&);
void new_nle3d01(unique_ptr<Material>&, istringstream&);
void new_orthotropicelastic3d(unique_ptr<Material>&, istringstream&);
void new_paraboliccc(unique_ptr<Material>&, istringstream&);
void new_polyelastic1d(unique_ptr<Material>&, istringstream&);
void new_polyj2(unique_ptr<Material>&, istringstream&);
void new_rambergosgood(unique_ptr<Material>&, istringstream&);
void new_rebar2d(unique_ptr<Material>&, istringstream&);
void new_rebar3d(unique_ptr<Material>&, istringstream&);
void new_sliplock(unique_ptr<Material>&, istringstream&);
void new_simplesand(unique_ptr<Material>&, istringstream&);
void new_steelbrb(unique_ptr<Material>&, istringstream&);
void new_tablecdp(unique_ptr<Material>&, istringstream&);
void new_tablegurson(unique_ptr<Material>&, istringstream&);
void new_trivial(unique_ptr<Material>&, istringstream&);
void new_vafcrp(unique_ptr<Material>&, istringstream&);
void new_vafcrp1d(unique_ptr<Material>&, istringstream&);
void new_viscosity01(unique_ptr<Material>&, istringstream&);
void new_viscosity02(unique_ptr<Material>&, istringstream&);
void new_bilinearviscosity(unique_ptr<Material>&, istringstream&);
void new_yeoh(unique_ptr<Material>&, istringstream&);

// wrapper
void new_axisymmetric(unique_ptr<Material>&, istringstream&);
void new_parallel(unique_ptr<Material>&, istringstream&);
void new_planestrain(unique_ptr<Material>&, istringstream&, unsigned);
void new_planestress(unique_ptr<Material>&, istringstream&);
void new_sequential(unique_ptr<Material>&, istringstream&);
void new_uniaxial(unique_ptr<Material>&, istringstream&);
void new_laminated(unique_ptr<Material>&, istringstream&);
void new_stacked(unique_ptr<Material>&, istringstream&);
void new_rotation2d(unique_ptr<Material>&, istringstream&);
void new_rotation3d(unique_ptr<Material>&, istringstream&);
void new_substepping(unique_ptr<Material>&, istringstream&);

// degradation
void new_trilineardegradation(unique_ptr<Material>&, istringstream&);
void new_dhakal(unique_ptr<Material>&, istringstream&);

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

#endif
