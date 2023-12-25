/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class Criterion
 * @brief A Criterion class.
 *
 * The Criterion class.
 *
 * @author tlc
 * @date 26/09/2017
 * @version 0.1.0
 * @file Criterion.h
 * @addtogroup Criterion
 * @{
 */

#ifndef CRITERION_H
#define CRITERION_H

#include <Domain/Tag.h>

class DomainBase;

class Criterion : public Tag {
    unsigned step_tag;

public:
    explicit Criterion(unsigned = 0, unsigned = 0);
    ~Criterion() override = default;

    void set_step_tag(unsigned);
    [[nodiscard]] unsigned get_step_tag() const;

    virtual unique_ptr<Criterion> get_copy() = 0;

    virtual int initialize(const shared_ptr<DomainBase>&);

    virtual int process(const shared_ptr<DomainBase>&) = 0;
};

#endif

//! @}
