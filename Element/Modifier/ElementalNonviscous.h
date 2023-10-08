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
 * @class ElementalNonviscous
 * @brief A ElementalNonviscous damping class.
 * @author tlc
 * @date 02/10/2023
 * @version 0.2.0
 * @file ElementalNonviscous.h
 * @addtogroup Modifier
 * @{
 */

#ifndef ELEMENTALNONVISCOUS_H
#define ELEMENTALNONVISCOUS_H

#include <Element/Modifier/Modifier.h>
#include <Domain/Factory.hpp>

class ElementalNonviscous : public Modifier {
    const cx_vec m, s;

    weak_ptr<Factory<double>> factory;

public:
    ElementalNonviscous(unsigned, cx_vec&&, cx_vec&&, uvec&& = {});

    [[nodiscard]] bool has_nonviscous() const override { return true; }

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;
};

class ElementalNonviscousGroup final : public ElementalNonviscous {
    const unsigned group_tag;

public:
    ElementalNonviscousGroup(unsigned, cx_vec&&, cx_vec&&, unsigned);

    int initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
