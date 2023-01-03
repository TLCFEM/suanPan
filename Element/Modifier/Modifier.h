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
 * @class Modifier
 * @brief A Modifier class.
 * @author tlc
 * @date 20/12/2018
 * @version 0.2.0
 * @file Modifier.h
 * @addtogroup Modifier
 * @{
 */

#ifndef MODIFIER_H
#define MODIFIER_H

#include <Domain/Tag.h>
#include <Element/Element.h>

class DomainBase;

class Modifier : public Tag {
    uvec element_tag;

protected:
    std::vector<weak_ptr<Element>> element_pool;

public:
    explicit Modifier(unsigned = 0, // tag
                      uvec&& = {}   // element tags
        );
    Modifier(const Modifier&) = delete;            // copy forbidden
    Modifier(Modifier&&) = delete;                 // move forbidden
    Modifier& operator=(const Modifier&) = delete; // assign forbidden
    Modifier& operator=(Modifier&&) = delete;      // assign forbidden

    ~Modifier() override = default;

    virtual void initialize(const shared_ptr<DomainBase>&);

    virtual int update_status() = 0;
};

#endif

//! @}
