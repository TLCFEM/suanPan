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
 * @class Expression
 * @brief A Expression class represents a maths expression.
 *
 * @author tlc
 * @date 12/01/2023
 * @file Expression.h
 * @addtogroup Utility
 * @{
 */

#ifndef EXPRESSION_H
#define EXPRESSION_H

#include <Domain/Tag.h>
#include <exprtk/exprtk.hpp>

class Expression final : public Tag {
    static std::mutex parser_lock;
    static exprtk::parser<double> parser;

    Col<double> x;

    exprtk::symbol_table<double> symbol_table;

    exprtk::expression<double> expression;

    std::string expression_text;

public:
    explicit Expression(unsigned, const std::string&);

    [[nodiscard]] uword size() const;

    bool compile(const std::string&);

    double evaluate(const Col<double>&);

    Col<double> gradient(const Col<double>&);

    void print() override;
};

#endif

//! @}
