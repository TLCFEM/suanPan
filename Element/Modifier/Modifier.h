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

class Modifier : public UniqueTag {
protected:
    uvec element_tag;

    std::vector<std::weak_ptr<Element>> element_pool;

public:
    explicit Modifier(
        unsigned = 0, // tag
        uvec&& = {}   // element tags
    );

    [[nodiscard]] virtual bool has_nonviscous() const { return false; }

    virtual int initialize(const shared_ptr<DomainBase>&);

    virtual bool if_apply(const shared_ptr<DomainBase>&);

    virtual int update_status() = 0;
};

class ModifierDynamics : public Modifier {
public:
    using Modifier::Modifier;

    bool if_apply(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
