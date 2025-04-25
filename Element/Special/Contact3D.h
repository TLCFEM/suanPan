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
 * @class Contact3D
 * @brief A Contact3D class.
 *
 * The Contact3D class defines node to triangular facet contact via penalty function.
 *
 * @author tlc
 * @date 29/04/2021
 * @version 0.1.0
 * @file Contact3D.h
 * @addtogroup Constraint
 * @{
 */

#ifndef CONTACT3D_H
#define CONTACT3D_H

#include <Element/Element.h>
#include <array>

class Contact3D final : public Element {
    struct SlaveNode {
        std::weak_ptr<Node> node;
        uvec local_span;
        vec position;
    };

    struct MasterNode {
        std::weak_ptr<Node> node;
        uvec local_span;
        vec position;
        vec outer_norm;
    };

    struct MasterFacet {
        std::array<MasterNode, 3> node;
        vec facet_outer_norm;
        double facet_area;
    };

    static constexpr unsigned c_dof = 3;

    const unsigned master_tag, slave_tag;

    std::vector<MasterFacet> master;
    std::vector<SlaveNode> slave;

    const double alpha;

    void update_position();
    void check_contact(const MasterFacet&, const SlaveNode&);

public:
    Contact3D(unsigned, unsigned, unsigned, double = 1E8);

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
