/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
 * @class Group
 * @brief The Group class.
 *
 * @author tlc
 * @date 21/03/2020
 * @version 0.1.0
 * @file Group.h
 * @addtogroup Group
 * @{
 */

#ifndef GROUP_H
#define GROUP_H

#include <Domain/Tag.h>

class DomainBase;

class Group : public Tag {
protected:
    uvec pool;

public:
    explicit Group(unsigned);
    Group(unsigned, uvec&&);
    Group(const Group&) = delete;            // copy forbidden
    Group(Group&&) = delete;                 // move forbidden
    Group& operator=(const Group&) = delete; // assign forbidden
    Group& operator=(Group&&) = delete;      // assign forbidden
    ~Group() override = default;

    virtual void initialize(const shared_ptr<DomainBase>&);

    [[nodiscard]] const uvec& get_pool() const;

    void print() override;
};

#endif

//! @}
