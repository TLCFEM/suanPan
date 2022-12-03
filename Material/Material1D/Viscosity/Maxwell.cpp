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

#include "Maxwell.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Recorder/OutputType.h>

Maxwell::Maxwell(const unsigned T, const unsigned DT, const unsigned ST, const bool UM, const unsigned PC, const double BT)
    : Material1D(T, 0.)
    , damper_tag(DT)
    , spring_tag(ST)
    , proceed(PC)
    , use_matrix(UM)
    , beta(std::min(std::max(0., std::fabs(BT)), 1.)) { access::rw(tolerance) = 1E-12; }

Maxwell::Maxwell(const Maxwell& old_obj)
    : Material1D(old_obj)
    , incre_time(old_obj.incre_time)
    , damper_tag(old_obj.damper_tag)
    , spring_tag(old_obj.spring_tag)
    , proceed(old_obj.proceed)
    , use_matrix(old_obj.use_matrix)
    , beta(old_obj.beta)
    , damper(suanpan::make_copy(old_obj.damper))
    , spring(suanpan::make_copy(old_obj.spring)) {}

int Maxwell::initialize(const shared_ptr<DomainBase>& D) {
    damper = suanpan::initialized_material_copy(D, damper_tag);
    spring = suanpan::initialized_material_copy(D, spring_tag);

    if(nullptr == damper || nullptr == spring) return SUANPAN_FAIL;

    incre_time = &get_incre_time(D->get_factory());

    trial_strain_rate = current_strain_rate = incre_strain_rate.zeros(1);

    trial_damping = current_damping = initial_damping = 0.;
    trial_stiffness = current_stiffness = initial_stiffness = 0.;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Maxwell::get_copy() { return make_unique<Maxwell>(*this); }

int Maxwell::update_trial_status(const vec&) {
    suanpan_error("Maxwell receives strain only from the associated element, check the model.\n");
    return SUANPAN_FAIL;
}

int Maxwell::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    incre_strain = (trial_strain = t_strain) - current_strain;
    incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

    if(fabs(incre_strain(0)) + fabs(incre_strain_rate(0)) <= datum::eps) return SUANPAN_SUCCESS;

    const auto& K1 = spring->get_trial_stiffness().at(0);
    const auto& K2 = damper->get_trial_stiffness().at(0);
    const auto& K3 = damper->get_trial_damping().at(0);
    const auto& F1 = spring->get_trial_stress().at(0);
    const auto& F2 = damper->get_trial_stress().at(0);

    // \beta\Delta{}t
    const auto factor_a = beta * *incre_time;

    const auto target = *incre_time * (current_strain_rate(0) - damper->get_current_strain_rate().at(0)) + factor_a * incre_strain_rate(0);

    vec solution(3, fill::zeros);

    counter = 0;

    if(double error, ref_error = 1.; use_matrix) {
        mat inv_jacobian(3, 3);

        inv_jacobian(0, 2) = -factor_a;
        inv_jacobian(1, 2) = factor_a;
        inv_jacobian(2, 2) = 1.;

        while(++counter < max_iteration) {
            const vec residual{incre_strain(0) - solution(0) - solution(1), target - solution(0) - factor_a * solution(2), F1 - F2};

            inv_jacobian(0, 0) = factor_a * K2;
            inv_jacobian(1, 0) = factor_a * K1 + K3;
            inv_jacobian(2, 0) = -K2;

            inv_jacobian(0, 1) = K3;
            inv_jacobian(1, 1) = -K3;
            inv_jacobian(2, 1) = K1 + K2;

            const vec incre = inv_jacobian * residual / (factor_a * (K1 + K2) + K3);

            if(1 == counter) ref_error = std::max(1., norm(residual));
            suanpan_extra_debug("Maxwell local iteration error: %.4E.\n", error = norm(residual) / ref_error);
            if(norm(incre) <= tolerance && error <= tolerance) break;
            solution += incre;
            spring->update_incre_status(solution(0));
            damper->update_incre_status(solution(1), solution(2));
        }
    }
    else
        while(++counter < max_iteration) {
            const auto residual_a = incre_strain(0) - solution(0) - solution(1);
            const auto residual_b = target - solution(0) - factor_a * solution(2);
            const auto residual_c = F1 - F2;
            const auto residual = residual_a * K2 - residual_c + residual_b / factor_a * K3;
            const auto jacobian = K1 + K2 + K3 / factor_a;
            const auto incre = residual / jacobian;
            if(1 == counter) ref_error = std::max(1., fabs(residual));
            suanpan_extra_debug("Maxwell local iteration error: %.4E.\n", error = fabs(residual) / ref_error);
            if(fabs(incre) <= tolerance && error <= tolerance) break;
            solution(0) += incre;
            solution(1) += residual_a - incre;
            solution(2) += (residual_b - incre) / factor_a;
            spring->update_incre_status(solution(0));
            damper->update_incre_status(solution(1), solution(2));
        }

    if(max_iteration != counter) {
        delay_counter = 0;

        trial_stress = .5 * (F1 + F2);

        trial_damping = trial_stiffness = factor_a / (factor_a * (K1 + K2) + K3) * K1;
        trial_stiffness *= K2;
        trial_damping *= K3;

        return SUANPAN_SUCCESS;
    }

    if(1 >= proceed || ++delay_counter == proceed) {
        suanpan_error("Maxwell local iteration cannot converge within %u iterations.\n", max_iteration);
        return SUANPAN_FAIL;
    }

    return reset_status();
}

int Maxwell::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_strain_rate = current_strain_rate.zeros();
    trial_damping = current_damping = initial_damping;
    trial_stiffness = current_stiffness = initial_stiffness;
    return spring->clear_status() + damper->clear_status();
}

int Maxwell::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_strain_rate = trial_strain_rate;
    current_damping = trial_damping;
    current_stiffness = trial_stiffness;
    return spring->commit_status() + damper->commit_status();
}

int Maxwell::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_strain_rate = current_strain_rate;
    trial_damping = current_damping;
    trial_stiffness = current_stiffness;
    return spring->reset_status() + damper->reset_status();
}

vector<vec> Maxwell::record(const OutputType P) {
    vector<vec> data;

    if(OutputType::SD == P || OutputType::SS == P || OutputType::S == P) data.emplace_back(current_stress);
    else if(OutputType::ED == P) data.emplace_back(damper->get_current_strain());
    else if(OutputType::VD == P) data.emplace_back(damper->get_current_strain_rate());
    else if(OutputType::ES == P) data.emplace_back(spring->get_current_strain());
    else if(OutputType::VS == P) data.emplace_back(current_strain_rate - damper->get_current_strain_rate());
    else if(OutputType::E == P) data.emplace_back(current_strain);
    else if(OutputType::V == P) data.emplace_back(current_strain_rate);
    else if(OutputType::LITR == P) data.emplace_back(vec{static_cast<double>(counter)});

    return data;
}

void Maxwell::print() { suanpan_info("A Maxwell material model.\n"); }
