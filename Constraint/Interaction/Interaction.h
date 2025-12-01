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
 * @class Interaction
 * @author tlc
 * @date 01/12/2025
 * @version 0.1.0
 * @file Interaction.h
 * @addtogroup Interaction
 * @{
 */

#ifndef INTERACTION_H
#define INTERACTION_H

#include <Domain/Tag.h>

class DomainBase;
class Element;

struct InteractionPair {
    shared_ptr<Element> object_i, object_j;

    InteractionPair(const shared_ptr<Element>&, const shared_ptr<Element>&);

    [[nodiscard]] double compression() const;
};

class Interaction : public CopyableTag {
    shared_ptr<DomainBase> domain;

public:
    Interaction() = default;

    void initialize(const shared_ptr<DomainBase>&);

    [[nodiscard]] virtual int apply(const shared_ptr<InteractionPair>&) const = 0;
};

#endif

//! @}
