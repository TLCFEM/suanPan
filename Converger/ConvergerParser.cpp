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

#include "ConvergerParser.h"

#include <Converger/Converger>
#include <Domain/DomainBase.h>
#include <Step/Step.h>
#include <Toolbox/utility.h>

int create_new_converger(const shared_ptr<DomainBase>& domain, istringstream& command) {
    std::string converger_id;
    if(!get_input(command, converger_id)) {
        suanpan_error("A valid converger type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    if(is_equal(converger_id.substr(0, 5), "Logic")) {
        unsigned tag_a, tag_b;
        if(!get_input(command, tag_a) || !get_input(command, tag_b)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(converger_id, "LogicAND") && domain->insert(std::make_shared<LogicAND>(tag, tag_a, tag_b))) code = 1; // NOLINT(bugprone-branch-clone)
        else if(is_equal(converger_id, "LogicOR") && domain->insert(std::make_shared<LogicOR>(tag, tag_a, tag_b))) code = 1;
        else if(is_equal(converger_id, "LogicXOR") && domain->insert(std::make_shared<LogicXOR>(tag, tag_a, tag_b))) code = 1;
        else
            suanpan_error("Cannot identify converger type.\n");
    }
    else {
        auto tolerance = 1E-6;
        if(!is_equal(converger_id, "FixedNumber") && (!command.eof() && !get_input(command, tolerance))) {
            suanpan_error("A valid tolerance is required.\n");
            return SUANPAN_SUCCESS;
        }

        auto max_iteration = 10;
        if(!command.eof() && !get_input(command, max_iteration)) {
            suanpan_error("A valid maximum iteration is required.\n");
            return SUANPAN_SUCCESS;
        }

        std::string print_flag = "false";
        if(!command.eof() && !get_input(command, print_flag)) {
            suanpan_error("A valid print flag is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(converger_id, "AbsResidual") && domain->insert(std::make_shared<AbsResidual>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1; // NOLINT(bugprone-branch-clone)
        else if(is_equal(converger_id, "RelResidual") && domain->insert(std::make_shared<RelResidual>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsIncreDisp") && domain->insert(std::make_shared<AbsIncreDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelIncreDisp") && domain->insert(std::make_shared<RelIncreDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsIncreAcc") && domain->insert(std::make_shared<AbsIncreAcc>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelIncreAcc") && domain->insert(std::make_shared<RelIncreAcc>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsDisp") && domain->insert(std::make_shared<AbsDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelDisp") && domain->insert(std::make_shared<RelDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsError") && domain->insert(std::make_shared<AbsError>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelError") && domain->insert(std::make_shared<RelError>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "AbsIncreEnergy") && domain->insert(std::make_shared<AbsIncreEnergy>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "RelIncreEnergy") && domain->insert(std::make_shared<RelIncreEnergy>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
        else if(is_equal(converger_id, "FixedNumber") && domain->insert(std::make_shared<FixedNumber>(tag, max_iteration, is_true(print_flag)))) code = 1;
        else
            suanpan_error("Cannot identify converger type.\n");
    }

    if(1 == code) {
        if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_converger_tag(tag);
        domain->set_current_converger_tag(tag);
    }
    else
        suanpan_error("Fail to create new converger via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
