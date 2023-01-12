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

#include "Expression.h"
#include <Toolbox/utility.h>

std::mutex Expression::parser_lock;
exprtk::parser<double> Expression::parser;

Expression::Expression(const unsigned tag, const std::string& variable_string)
    : Tag(tag) {
    const auto variable_list = suanpan::expression::split(variable_string);

    x.zeros(variable_list.size());
    for(size_t I = 0; I < variable_list.size(); ++I) symbol_table.add_variable(variable_list[I], x(I));

    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);
}

uword Expression::size() const { return x.n_elem; }

bool Expression::compile(const std::string& expression_string) {
    expression_text = expression_string;
    const auto lock = std::scoped_lock{parser_lock};
    return parser.compile(expression_string, expression);
}

double Expression::evaluate(const Col<double>& in_x) {
    x = in_x;
    return expression.value();
}

Col<double> Expression::gradient(const Col<double>& in_x) {
    x = in_x;
    Col<double> result(x.n_elem);
    for(uword I = 0; I < x.n_elem; ++I) result(I) = derivative(expression, x(I));
    return result;
}

void Expression::print() {
    suanpan_info("An expression represents \"{}\".", expression_text);
}
