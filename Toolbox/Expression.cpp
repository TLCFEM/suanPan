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

std::mutex Expression::parser_mutex;
exprtk::parser<double> Expression::parser;

Expression::Expression(const unsigned tag, const std::string& variable_string)
    : Tag(tag) {
    const auto variable_list = suanpan::expression::split(variable_string);

    std::vector<std::pair<string, unsigned>> variable_size_list;
    for(const auto& I : variable_list)
        if(is_integer(I)) variable_size_list.back().second = std::stoi(I);
        else variable_size_list.emplace_back(I, 1);

    x.zeros(std::accumulate(variable_size_list.cbegin(), variable_size_list.cend(), 0, [](const auto a, const auto b) { return a + b.second; }));

    unsigned index = 0;
    for(const auto& [name, size] : variable_size_list)
        if(1 == size) symbol_table.add_variable(name, x(index++));
        else {
            symbol_table.add_vector(name, x.memptr() + index, size);
            index += size;
        }

    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);
}

uword Expression::size() const { return x.n_elem; }

bool Expression::compile(const std::string& expression_string) {
    expression_text = expression_string;
    std::scoped_lock lock{parser_mutex};
    return parser.compile(expression_string, expression);
}

string Expression::error() { return parser.error(); }

double Expression::evaluate(const double in_x) { return evaluate(vec{in_x}); }

Col<double> Expression::gradient(const double in_x) { return gradient(vec{in_x}); }

void Expression::print() {
    suanpan_info("An expression represents \"{}\".", expression_text);
}

double SimpleScalarExpression::evaluate(const Col<double>& in_x) {
    suanpan_assert([&] { if(x.n_elem != in_x.n_elem) throw std::runtime_error("input size mismatch"); });

    x = in_x;
    return expression.value();
}

Col<double> SimpleScalarExpression::gradient(const Col<double>& in_x) {
    suanpan_assert([&] { if(x.n_elem != in_x.n_elem) throw std::runtime_error("input size mismatch"); });

    x = in_x;
    Col<double> result(x.n_elem);
    for(uword I = 0; I < x.n_elem; ++I) result(I) = derivative(expression, x(I));
    return result;
}
