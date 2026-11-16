/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

namespace {
    template<typename T> void new_logic_converger(unique_ptr<Converger>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned tag_a, tag_b;
        if(!get_input(command, tag_a, tag_b)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag, tag_a, tag_b);
    }

    template<typename T> void new_default_converger(unique_ptr<Converger>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto tolerance = 1E-6;
        if(!command.eof() && !get_input(command, tolerance)) {
            suanpan_error("A valid tolerance is required.\n");
            return;
        }

        auto max_iteration = 10;
        if(!command.eof() && !get_input(command, max_iteration)) {
            suanpan_error("A valid maximum iteration is required.\n");
            return;
        }

        std::string print_flag = "false";
        if(!command.eof() && !get_input(command, print_flag)) {
            suanpan_error("A valid print flag is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag, tolerance, max_iteration, is_true(print_flag));
    }

    void new_fixed_number(unique_ptr<Converger>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto max_iteration = 10;
        if(!command.eof() && !get_input(command, max_iteration)) {
            suanpan_error("A valid maximum iteration is required.\n");
            return;
        }

        std::string print_flag = "false";
        if(!command.eof() && !get_input(command, print_flag)) {
            suanpan_error("A valid print flag is required.\n");
            return;
        }

        return_obj = std::make_unique<FixedNumber>(tag, max_iteration, is_true(print_flag));
    }
} // namespace

int create_new_converger(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string converger_id;
    if(!get_input(command, converger_id)) {
        suanpan_error("A valid converger type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Converger> new_converger = nullptr;

    if(is_equal_any(converger_id, "LogicAND", "LogicConvergerAND")) new_logic_converger<LogicConvergerAND>(new_converger, command);
    else if(is_equal_any(converger_id, "LogicOR", "LogicConvergerOR")) new_logic_converger<LogicConvergerOR>(new_converger, command);
    else if(is_equal_any(converger_id, "LogicXOR", "LogicConvergerXOR")) new_logic_converger<LogicConvergerXOR>(new_converger, command);
    else if(is_equal(converger_id, "AbsResidual")) new_default_converger<AbsResidual>(new_converger, command);
    else if(is_equal(converger_id, "RelResidual")) new_default_converger<RelResidual>(new_converger, command);
    else if(is_equal(converger_id, "AbsIncreDisp")) new_default_converger<AbsIncreDisp>(new_converger, command);
    else if(is_equal(converger_id, "RelIncreDisp")) new_default_converger<RelIncreDisp>(new_converger, command);
    else if(is_equal(converger_id, "AbsIncreAcc")) new_default_converger<AbsIncreAcc>(new_converger, command);
    else if(is_equal(converger_id, "RelIncreAcc")) new_default_converger<RelIncreAcc>(new_converger, command);
    else if(is_equal(converger_id, "AbsDisp")) new_default_converger<AbsDisp>(new_converger, command);
    else if(is_equal(converger_id, "RelDisp")) new_default_converger<RelDisp>(new_converger, command);
    else if(is_equal(converger_id, "AbsError")) new_default_converger<AbsError>(new_converger, command);
    else if(is_equal(converger_id, "RelError")) new_default_converger<RelError>(new_converger, command);
    else if(is_equal(converger_id, "AbsIncreEnergy")) new_default_converger<AbsIncreEnergy>(new_converger, command);
    else if(is_equal(converger_id, "RelIncreEnergy")) new_default_converger<RelIncreEnergy>(new_converger, command);
    else if(is_equal(converger_id, "FixedNumber")) new_fixed_number(new_converger, command);

    if(nullptr == new_converger) {
        suanpan_error("Fail to create new converger via \"{}\".\n", command.str());
        return SUANPAN_SUCCESS;
    }

    const auto tag = new_converger->get_tag();

    if(!domain->insert(std::move(new_converger))) {
        suanpan_error("Fail to create new converger via \"{}\".\n", command.str());
        return SUANPAN_SUCCESS;
    }

    if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_converger_tag(tag);
    domain->set_current_converger_tag(tag);

    return SUANPAN_SUCCESS;
}
