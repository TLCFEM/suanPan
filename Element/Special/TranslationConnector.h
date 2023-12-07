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
 * @class TranslationConnector
 * @brief A TranslationConnector class.
 *
 * The TranslationConnector class.
 *
 * @author tlc
 * @date 06/03/2023
 * @version 0.1.0
 * @file TranslationConnector.h
 * @addtogroup Constraint
 * @{
 */

#ifndef TRANSLATIONCONNECTOR_H
#define TRANSLATIONCONNECTOR_H

#include <Element/Element.h>

class TranslationConnector : public Element {
    static const unsigned c_node = 3;

    const unsigned c_dof;

    const unsigned c_size = c_node * c_dof;

    const span sa, sb, sc;

    const double alpha;

    const double s_near{0.}, s_far{0.};

public:
    TranslationConnector(unsigned, uvec&&, unsigned, double);

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int clear_status() override { return SUANPAN_SUCCESS; }

    int commit_status() override { return SUANPAN_SUCCESS; }

    int reset_status() override { return SUANPAN_SUCCESS; }
};

class TranslationConnector2D final : public TranslationConnector {
public:
    TranslationConnector2D(const unsigned T, uvec&& N, const double P)
        : TranslationConnector(T, std::move(N), 2u, P) {}
};

class TranslationConnector3D final : public TranslationConnector {
public:
    TranslationConnector3D(const unsigned T, uvec&& N, const double P)
        : TranslationConnector(T, std::move(N), 3u, P) {}
};

#endif

//! @}
