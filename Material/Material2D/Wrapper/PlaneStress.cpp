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

#include "PlaneStress.h"
#include <Domain/DomainBase.h>

const uvec PlaneStress::F1{0, 1, 3};
const uvec PlaneStress::F2{2, 4, 5};

mat PlaneStress::form_stiffness(const mat& full_stiffness) {
    const auto& A = full_stiffness(F1, F1);
    const auto& U = full_stiffness(F1, F2);
    const auto& C = full_stiffness(F2, F2);
    const auto& V = full_stiffness(F2, F1);

    return A - U * solve(C, V);

    /*
    const auto S = inv(A);
    const auto VS = V * S;

    return inv(S * (eye(3, 3) - U * solve(VS * U - C, VS)));
    */
}

PlaneStress::PlaneStress(const unsigned T, const unsigned BT, const unsigned MI, const bool FM)
    : Material2D(T, PlaneType::S, 0.)
    , base_tag(BT)
    , max_iteration(MI)
    , use_full_matrix(FM) { access::rw(tolerance) = 1E-12; }

PlaneStress::PlaneStress(const PlaneStress& old_obj)
    : Material2D(old_obj)
    , base_tag(old_obj.base_tag)
    , max_iteration(old_obj.max_iteration)
    , use_full_matrix(old_obj.use_full_matrix)
    , base(suanpan::make_copy(old_obj.base))
    , trial_full_strain(old_obj.trial_full_strain)
    , current_full_strain(old_obj.current_full_strain) {}

int PlaneStress::initialize(const shared_ptr<DomainBase>& D) {
    base = suanpan::initialized_material_copy(D, base_tag);

    if(nullptr == base || base->get_material_type() != MaterialType::D3) {
        SP_E("A valid 3D host material is required.\n");
        return SUANPAN_FAIL;
    }

    trial_full_strain = current_full_strain.zeros(6);

    current_stiffness = trial_stiffness = initial_stiffness = form_stiffness(base->get_initial_stiffness());

    return SUANPAN_SUCCESS;
}

double PlaneStress::get_parameter(const ParameterType P) const {
    if(ParameterType::PLANETYPE == P) return static_cast<double>(plane_type);
    return base->get_parameter(P);
}

unique_ptr<Material> PlaneStress::get_copy() { return make_unique<PlaneStress>(*this); }

int PlaneStress::update_trial_status(const vec& t_strain) {
    auto& t_stress = base->get_trial_stress();
    auto& t_stiffness = base->get_trial_stiffness();

    if(norm(incre_strain = t_strain - trial_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_full_strain(F1) = trial_strain = t_strain;

    if(1 != max_iteration) {
        unsigned counter = 0;

        if(use_full_matrix)
            while(++counter < max_iteration) {
                if(SUANPAN_SUCCESS != base->update_trial_status(trial_full_strain)) return SUANPAN_FAIL;
                trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2));
                const auto error = norm(t_stress(F2));
                SP_D("Local iteration error: {:.5E}.\n", error);
                if(error < tolerance) break;
            }
        else
            while(++counter < max_iteration) {
                if(SUANPAN_SUCCESS != base->update_trial_status(trial_full_strain)) return SUANPAN_FAIL;
                trial_full_strain(F2) -= t_stress(F2) / vec(t_stiffness.diag())(F2);
                const auto error = norm(t_stress(F2));
                SP_D("Local iteration error: {:.5E}.\n", error);
                if(error < tolerance) break;
            }

        if(counter == max_iteration) {
            SP_E("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        trial_stress = t_stress(F1);
    }
    else {
        trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2) + t_stiffness(F2, F1) * incre_strain);

        if(SUANPAN_SUCCESS != base->update_trial_status(trial_full_strain)) return SUANPAN_FAIL;

        trial_stress = t_stress(F1) - t_stiffness(F1, F2) * solve(t_stiffness(F2, F2), t_stress(F2));
    }

    trial_stiffness = form_stiffness(t_stiffness);

    return SUANPAN_SUCCESS;
}

int PlaneStress::clear_status() {
    current_full_strain.zeros();
    trial_full_strain.zeros();
    current_strain.zeros();
    trial_strain.zeros();
    current_stress.zeros();
    trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return base->clear_status();
}

int PlaneStress::commit_status() {
    current_full_strain = trial_full_strain;
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return base->commit_status();
}

int PlaneStress::reset_status() {
    trial_full_strain = current_full_strain;
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return base->reset_status();
}

vector<vec> PlaneStress::record(const OutputType P) { return base->record(P); }

void PlaneStress::print() {
    sp_info("A plane stress wrapper.\n");
    sp_info("Strain:", current_strain);
    sp_info("Stress:", current_stress);
    if(base) base->print();
}
