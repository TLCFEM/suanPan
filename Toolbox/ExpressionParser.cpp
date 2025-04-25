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

#include "ExpressionParser.h"
#include <Domain/DomainBase.h>
#include <Toolbox/Expression.h>
#include <Toolbox/utility.h>

int check_file(std::string& expression) {
    if(!fs::exists(expression)) return SUANPAN_SUCCESS;

    const std::ifstream file(expression);
    if(!file.is_open()) {
        suanpan_error("Fail to open \"{}\".\n", expression);
        return SUANPAN_FAIL;
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    expression = buffer.str();

    return SUANPAN_SUCCESS;
}

void new_simplescalar(unique_ptr<Expression>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid expression tag is required.\n");
        return;
    }

    std::string variable_list;
    if(!get_input(command, variable_list)) {
        suanpan_error("A valid vertical bar separated variable list is required.\n");
        return;
    }

    std::string expression;
    if(!get_input(command, expression)) {
        suanpan_error("A valid expression or expression file name is required.\n");
        return;
    }

    if(SUANPAN_SUCCESS != check_file(expression)) return;

    return_obj = make_unique<SimpleScalarExpression>(tag, variable_list);

    if(!return_obj->compile(expression)) {
        suanpan_error("Fail to parse \"{}\", error: {}.\n", expression, return_obj->error());
        return_obj = nullptr;
    }
}

void new_simplevector(unique_ptr<Expression>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid expression tag is required.\n");
        return;
    }

    std::string input_list, output_list;
    if(!get_input(command, input_list, output_list)) {
        suanpan_error("A valid vertical bar separated variable list is required.\n");
        return;
    }

    std::string expression;
    if(!get_input(command, expression)) {
        suanpan_error("A valid expression or expression file name is required.\n");
        return;
    }

    if(SUANPAN_SUCCESS != check_file(expression)) return;

    return_obj = make_unique<SimpleVectorExpression>(tag, input_list, output_list);

    if(!return_obj->compile(expression)) {
        suanpan_error("Fail to parse \"{}\", error: {}.\n", expression, return_obj->error());
        return_obj = nullptr;
    }
}

int create_new_expression(const shared_ptr<DomainBase>& domain, istringstream& command) {
    std::string expression_type;
    if(!get_input(command, expression_type)) {
        suanpan_error("A valid expression type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Expression> new_expression = nullptr;

    if(is_equal(expression_type, "SimpleScalar")) new_simplescalar(new_expression, command);
    else if(is_equal(expression_type, "SimpleVector")) new_simplevector(new_expression, command);

    if(nullptr == new_expression || !domain->insert(std::move(new_expression)))
        suanpan_error("Fail to create new expression via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
