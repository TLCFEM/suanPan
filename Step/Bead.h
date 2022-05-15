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
/**
 * @class Bead
 * @brief A Bead class is a top level container.
 * @author tlc
 * @date 01/10/2017
 * @version 0.3.0
 * @file Bead.h
 */

#ifndef BEAD_H
#define BEAD_H

#include <Domain/Storage.hpp>

class DomainBase;

class Bead {
    unsigned current_domain_tag = 1;

    DomainBaseStorage domain_pool;

public:
    Bead();

    bool insert(const shared_ptr<DomainBase>&);
    void erase_domain(unsigned);
    void enable_domain(unsigned);
    void disable_domain(unsigned);

    void set_current_domain_tag(unsigned);
    unsigned get_current_domain_tag() const;

    const shared_ptr<DomainBase>& get_domain(unsigned) const;
    const shared_ptr<DomainBase>& get_current_domain() const;

    friend shared_ptr<DomainBase>& get_domain(const shared_ptr<Bead>&, unsigned);
    friend shared_ptr<DomainBase>& get_current_domain(const shared_ptr<Bead>&);

    int precheck();

    int analyze();
};

#endif
