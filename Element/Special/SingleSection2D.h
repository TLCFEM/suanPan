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
 * @class SingleSection2D
 * @brief The SingleSection2D class.
 *
 * @author tlc
 * @date 06/10/2017
 * @version 0.1.0
 * @file SingleSection2D.h
 * @addtogroup Special
 * @ingroup Element
 * @{
 */

#ifndef SINGLESECTION2D_H
#define SINGLESECTION2D_H

#include <Element/SectionElement.h>

class SingleSection2D final : public SectionElement2D {
    static constexpr unsigned s_node = 1, s_dof = 2;

    unique_ptr<Section> s_section;

public:
    SingleSection2D(
        unsigned, // tag
        unsigned, // node tag
        unsigned  // section tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
