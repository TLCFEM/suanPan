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

#ifndef ELEMENTEXAMPLE_H
#define ELEMENTEXAMPLE_H

#include <Element/Element.h>

class ElementExample final : public Element {
    static constexpr unsigned m_node = 3, m_dof = 2, m_size = m_node * m_dof;

    const double thickness, area = 0.;

    mat strain_mat;

    unique_ptr<Material> m_material;

public:
    ElementExample(unsigned, uvec&&, unsigned, double = 1.);

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    void print() override;
};

#endif
