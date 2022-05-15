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
 * @class Embedded
 * @brief A Embedded class.
 *
 * The Embedded class.
 *
 * @author tlc
 * @date 12/08/2020
 * @version 0.1.0
 * @file Embedded.h
 * @addtogroup Constraint
 * @{
 */

#ifndef EMBEDDED_H
#define EMBEDDED_H

#include <Element/Element.h>

class Embedded : public Element {
    static constexpr unsigned max_iteration = 20;

    const unsigned e_dof;
    const unsigned host_tag;
    const unsigned host_size = 0;
    const double alpha;
    const rowvec iso_n;

    std::vector<uvec> idx;

    shared_ptr<Element> host_element;

public:
    Embedded(unsigned, unsigned, unsigned, unsigned, double);

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

class Embedded2D final : public Embedded {
public:
    Embedded2D(unsigned, unsigned, unsigned, double = 1E6);
};

class Embedded3D final : public Embedded {
public:
    Embedded3D(unsigned, unsigned, unsigned, double = 1E6);
};

#endif

//! @}
