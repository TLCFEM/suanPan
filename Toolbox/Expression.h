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

protected:
    Col<double> x;

    std::string expression_text;

    std::vector<std::string> variable_text_list;

    exprtk::expression<double> expression;

    exprtk::symbol_table<double> symbol_table;

public:
    Expression(unsigned, std::vector<std::string>&&);
    Expression(const Expression&) = delete;
    Expression(Expression&&) noexcept = delete;
    Expression& operator=(const Expression&) = delete;
    Expression& operator=(Expression&&) noexcept = delete;
    ~Expression() override = default;

    [[nodiscard]] virtual unique_ptr<Expression> get_copy() const = 0;

    [[nodiscard]] virtual uword input_size() const;
    [[nodiscard]] virtual uword output_size() const;

    bool compile(const std::string_view&);
    static std::string error();

    Mat<double> evaluate(double);
    virtual Mat<double> evaluate(const Col<double>&) = 0;

    Mat<double> gradient(double);
    virtual Mat<double> gradient(const Col<double>&) = 0;

    void print() override;
};

class SimpleScalarExpression final : public Expression {
public:
    SimpleScalarExpression(unsigned, const std::string_view&);

    [[nodiscard]] unique_ptr<Expression> get_copy() const override;

    Mat<double> evaluate(const Col<double>&) override;

    Mat<double> gradient(const Col<double>&) override;
};

class SimpleVectorExpression final : public Expression {
    Col<double> y;

public:
    SimpleVectorExpression(unsigned, const std::string_view&, const std::string_view&);

    [[nodiscard]] unique_ptr<Expression> get_copy() const override;

    [[nodiscard]] uword output_size() const override;

    Mat<double> evaluate(const Col<double>&) override;

    Mat<double> gradient(const Col<double>&) override { throw std::runtime_error("gradient is not implemented for vector expression"); }
};

#endif

//! @}
