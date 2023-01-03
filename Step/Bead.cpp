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

#include "Bead.h"
#include <Domain/Domain.h>
#include <Step/Step.h>

Bead::Bead() { insert(make_shared<Domain>(1)); }

bool Bead::insert(const shared_ptr<DomainBase>& D) { return domain_pool.insert(D); }

void Bead::erase_domain(const unsigned T) {
    if(domain_pool.find(T)) {
        domain_pool.erase(T);
        if(domain_pool.size() == 0) {
            insert(make_shared<Domain>(1));
            set_current_domain_tag(1);
            sp_info("Remove the last domain and switch to default one.\n");
        }
        else if(current_domain_tag == T) {
            set_current_domain_tag(domain_pool.cbegin()->first);
            sp_info("Switch to domain {}.\n", current_domain_tag);
        }
    }
    else sp_info("Cannot find domain {}, nothing changed.\n", T);
}

void Bead::enable_domain(const unsigned T) { domain_pool.enable(T); }

void Bead::disable_domain(const unsigned T) { domain_pool.disable(T); }

void Bead::set_current_domain_tag(const unsigned T) { current_domain_tag = T; }

unsigned Bead::get_current_domain_tag() const { return current_domain_tag; }

const shared_ptr<DomainBase>& Bead::get_domain(const unsigned T) const { return domain_pool.at(T); }

const shared_ptr<DomainBase>& Bead::get_current_domain() const { return domain_pool.at(current_domain_tag); }

int Bead::precheck() {
    for(const auto& [d_tag, t_domain] : domain_pool)
        if(t_domain->is_active())
            for(const auto& [s_tag, t_step] : t_domain->get_step_pool()) {
                t_domain->set_current_step_tag(t_step->get_tag());
                t_step->set_domain(t_domain);
                if(SUANPAN_SUCCESS != t_step->Step::initialize()) return SUANPAN_FAIL;
                if(SUANPAN_SUCCESS != t_step->initialize()) return SUANPAN_FAIL;
            }

    for(const auto& [d_tag, t_domain] : domain_pool) if(t_domain->is_active()) t_domain->clear_status();

    return SUANPAN_SUCCESS;
}

int Bead::analyze() {
    for(const auto& [d_tag, t_domain] : domain_pool) {
        if(!t_domain->is_active()) continue;
        bool initial_record = true;
        for(const auto& [s_tag, t_step] : t_domain->get_step_pool()) {
            t_domain->set_current_step_tag(t_step->get_tag());
            t_step->set_domain(t_domain);
            if(SUANPAN_FAIL == t_step->Step::initialize()) return SUANPAN_FAIL;
            if(SUANPAN_FAIL == t_step->initialize()) return SUANPAN_FAIL;
            if(initial_record) {
                initial_record = false;
                t_domain->record();
            }
            if(SUANPAN_FAIL == t_step->analyze()) return SUANPAN_FAIL;
        }
    }

    return SUANPAN_SUCCESS;
}

shared_ptr<DomainBase>& get_domain(const shared_ptr<Bead>& B, const unsigned T) { return B->domain_pool[T]; }

shared_ptr<DomainBase>& get_current_domain(const shared_ptr<Bead>& B) { return B->domain_pool[B->current_domain_tag]; }
