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

class Expression : public Tag {
    static std::mutex parser_mutex;
    static exprtk::parser<double> parser;

    std::string expression_text;

protected:
    Col<double> x;

    exprtk::expression<double> expression;

    exprtk::symbol_table<double> symbol_table;

public:
    explicit Expression(unsigned, const std::string&);

    [[nodiscard]] uword size() const;

    bool compile(const std::string&);
    static string error();

    Mat<double> evaluate(double);
    virtual Mat<double> evaluate(const Col<double>&) = 0;

    Mat<double> gradient(double);
    virtual Mat<double> gradient(const Col<double>&) = 0;

    void print() override;
};

class SimpleScalarExpression : public Expression {
public:
    using Expression::Expression;

    Mat<double> evaluate(const Col<double>&) override;

    Mat<double> gradient(const Col<double>&) override;
};

class SimpleVectorExpression : public Expression {
    Col<double> y;

public:
    SimpleVectorExpression(unsigned, const std::string&, const std::string&);

    Mat<double> evaluate(const Col<double>&) override;

    Mat<double> gradient(const Col<double>&) override { throw std::runtime_error("gradient is not implemented for vector expression"); }
};

#endif

//! @}
