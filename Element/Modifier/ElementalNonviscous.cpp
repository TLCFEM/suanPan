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

#include "ElementalNonviscous.h"
#include <Domain/DomainBase.h>
#include <Domain/Group/Group.h>

ElementalNonviscous::ElementalNonviscous(const unsigned T, cx_vec&& M, cx_vec&& S, uvec&& ET)
    : ModifierDynamics(T, std::move(ET))
      , m(std::move(M))
      , s(std::move(S)) {}

int ElementalNonviscous::initialize(const shared_ptr<DomainBase>& D) {
    ModifierDynamics::initialize(D);

    factory = D->get_factory();

    if(std::any_of(element_pool.cbegin(), element_pool.cend(), [](const weak_ptr<Element>& ele_ptr) { return !ele_ptr.lock()->get_current_nonviscous_force().empty(); })) {
        suanpan_error("Repeated element tags are detected, modifier {} is disabled.\n", get_tag());
        element_pool.clear();
        return SUANPAN_FAIL;
    }

    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        const auto t_ele = ele_ptr.lock();
        access::rw(t_ele->get_current_nonviscous_force()).zeros(t_ele->get_total_number(), m.n_elem);
    });

    return SUANPAN_SUCCESS;
}

int ElementalNonviscous::update_status() {
    // shortcircuit as computing complex parameters is expensive
    if(element_pool.empty()) return SUANPAN_SUCCESS;

    const cx_vec t_para = 2. / factory.lock()->get_incre_time() + s;
    const cx_vec s_para = (t_para - 2. * s) / t_para;
    const cx_vec m_para = m / t_para;
    const auto accu_para = accu(m_para).real();

    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& t_ptr) {
        const auto t_element = t_ptr.lock();

        if(t_element->get_current_nonviscous_force().n_cols != s_para.n_elem) return;

        auto& trial_nonviscous_force = access::rw(t_element->get_trial_nonviscous_force());
        trial_nonviscous_force = t_element->get_current_nonviscous_force() * diagmat(s_para) + (t_element->get_current_velocity() + t_element->get_trial_velocity()) * m_para.t();

        auto& trial_nonviscous = access::rw(t_element->get_trial_nonviscous());
        trial_nonviscous.zeros(t_element->get_total_number(), t_element->get_total_number());
        trial_nonviscous.diag().fill(accu_para);
    });

    return SUANPAN_SUCCESS;
}

ElementalNonviscousGroup::ElementalNonviscousGroup(const unsigned T, cx_vec&& M, cx_vec&& S, const unsigned GT)
    : ElementalNonviscous(T, std::move(M), std::move(S))
    , group_tag(GT) {}

int ElementalNonviscousGroup::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_group(group_tag)) return SUANPAN_FAIL;

    element_tag = unique(D->get_group(group_tag)->get_pool());

    return ElementalNonviscous::initialize(D);
}
