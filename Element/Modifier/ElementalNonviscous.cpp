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
#include <Domain/Factory.hpp>

ElementalNonviscous::ElementalNonviscous(const unsigned T, cx_vec&& M, cx_vec&& S, uvec&& ET)
    : Modifier(T, std::forward<uvec>(ET))
    , m(std::forward<cx_vec>(M))
    , s(std::forward<cx_vec>(S)) {}

int ElementalNonviscous::initialize(const shared_ptr<DomainBase>& D) {
    Modifier::initialize(D);

    if(std::ranges::any_of(element_pool, [](const weak_ptr<Element>& ele_ptr) { return !ele_ptr.lock()->get_current_nonviscous_force().empty(); })) {
        suanpan_error("Repeated element tags are detected.\n");
        return SUANPAN_FAIL;
    }

    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        const auto t_ele = ele_ptr.lock();
        access::rw(t_ele->get_current_nonviscous_force()).zeros(t_ele->get_total_number(), m.n_elem);
    });

    incre_time = &D->get_factory()->modify_incre_time();
    return SUANPAN_SUCCESS;
}

int ElementalNonviscous::update_status() {
    const cx_vec t_para = 2. / *incre_time + s;
    const cx_vec s_para = (t_para - 2. * s) / t_para;
    const cx_vec m_para = m / t_para;
    const auto accu_para = accu(m_para).real();

    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        const auto t_ele = ele_ptr.lock();

        auto& trial_nonviscous_force = access::rw(t_ele->get_trial_nonviscous_force());
        trial_nonviscous_force = t_ele->get_current_nonviscous_force() * diagmat(s_para) + (t_ele->get_current_velocity() + t_ele->get_trial_velocity()) * m_para.t();

        auto& trial_nonviscous = access::rw(t_ele->get_trial_nonviscous());
        trial_nonviscous.zeros(t_ele->get_total_number(), t_ele->get_total_number());
        trial_nonviscous.diag().fill(accu_para);
    });

    return SUANPAN_SUCCESS;
}
