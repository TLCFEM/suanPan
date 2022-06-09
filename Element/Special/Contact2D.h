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
 * @class Contact2D
 * @brief A Contact2D class.
 *
 * The Contact2D class.
 *
 * @author tlc
 * @date 29/06/2020
 * @version 0.1.0
 * @file Contact2D.h
 * @addtogroup Constraint
 * @{
 */

#ifndef CONTACT2D_H
#define CONTACT2D_H

#include <Element/Element.h>

class Contact2D final : public Element {
    struct MasterNode {
        weak_ptr<Node> node;
        vec position;
        vec axis;
        vec norm;
    };

    struct SlaveNode {
        weak_ptr<Node> node;
        vec position;
    };

    static constexpr unsigned c_dof = 2;

    static const mat rotation;

    const unsigned master_tag, slave_tag;

    std::vector<MasterNode> master;
    std::vector<SlaveNode> slave;

    const double alpha;

    void update_position();

public:
    Contact2D(unsigned, unsigned, unsigned, double = 1E8);

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
