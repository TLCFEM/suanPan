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

#include "Nonviscous01.h"
#include <Recorder/OutputType.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Step/Step.h>

Nonviscous01::Nonviscous01(const unsigned T, cx_vec&& M, cx_vec&& S)
    : DataNonviscous01{std::forward<cx_vec>(M), std::forward<cx_vec>(S)}
    , Material1D(T, 0.) {}

int Nonviscous01::initialize(const shared_ptr<DomainBase>& D) {
    if(nullptr != D) incre_time = &D->get_factory()->modify_incre_time();

    complex_damping.zeros(m.n_elem);
    s_para.zeros(m.n_elem);
    m_para.zeros(m.n_elem);

    const auto ini_time = D->get_current_step()->get_ini_step_size();
    trial_damping = current_damping = initial_damping = accu(ini_time * m / (2. + ini_time * s)).real();

    trial_strain_rate = current_strain_rate.zeros(1);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Nonviscous01::get_copy() { return make_unique<Nonviscous01>(*this); }

int Nonviscous01::update_trial_status(const vec&) {
    suanpan_error("Receives strain only from the associated element.\n");
    return SUANPAN_FAIL;
}

int Nonviscous01::update_trial_status(const vec&, const vec& t_strain_rate) {
    incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

    if(fabs(incre_strain_rate(0)) <= datum::eps) return SUANPAN_SUCCESS;

    const cx_vec t_para = 2. + *incre_time * s;
    s_para = (4. - t_para) / t_para;
    m_para = *incre_time * m / t_para;
    accu_para = accu(m_para).real();

    trial_stress = real(dot(complex_damping, s_para) + accu_para * (current_strain_rate + trial_strain_rate));

    trial_damping = accu_para;

    return SUANPAN_SUCCESS;
}

int Nonviscous01::clear_status() {
    complex_damping.zeros();
    s_para.zeros();
    m_para.zeros();
    current_strain_rate.zeros();
    current_stress.zeros();
    current_damping = initial_damping;
    return reset_status();
}

int Nonviscous01::commit_status() {
    complex_damping %= s_para;
    complex_damping += m_para * (trial_strain_rate + current_strain_rate);

    current_strain_rate = trial_strain_rate;
    current_stress = trial_stress;
    current_damping = trial_damping;
    return SUANPAN_SUCCESS;
}

int Nonviscous01::reset_status() {
    trial_strain_rate = current_strain_rate;
    trial_stress = current_stress;
    trial_damping = current_damping;
    return SUANPAN_SUCCESS;
}

vector<vec> Nonviscous01::record(const OutputType P) {
    if(OutputType::V == P) return {current_strain_rate};

    return Material1D::record(P);
}

void Nonviscous01::print() {
    suanpan_info("A uniaxial nonviscous damping material model.\n");
    Material1D::print();
}
