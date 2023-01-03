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
/**
 * @class Database
 * @brief A Database class is a top level container.
 * @author tlc
 * @date 27/08/2017
 * @version 0.2.1
 * @file Database.h
 */

#ifndef DATABASE_H
#define DATABASE_H

#include <Domain/Tag.h>

class DomainBase;

class Database : public Tag {
    shared_ptr<DomainBase> domain = nullptr;

public:
    explicit Database(unsigned = 0);
    ~Database() override = default;

    void set_domain(const shared_ptr<DomainBase>& D);
    [[nodiscard]] const shared_ptr<DomainBase>& get_domain() const;

    virtual int save() = 0;
};

#endif
