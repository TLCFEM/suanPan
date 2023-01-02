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

#include "LogicConverger.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <future>

LogicConverger::LogicConverger(const unsigned T, const unsigned TA, const unsigned TB)
    : Converger(T)
    , tag_a(TA)
    , tag_b(TB) {}

int LogicConverger::initialize() {
    const auto t_domain = get_domain().lock();
    if(nullptr == t_domain) {
        SP_E("A valid domain is required.\n");
        return SUANPAN_FAIL;
    }

    const auto& t_converger_a = t_domain->get_converger(tag_a);
    const auto& t_converger_b = t_domain->get_converger(tag_b);

    if(nullptr == t_converger_a || nullptr == t_converger_b) {
        SP_E("Two valid convergers are required.\n");
        return SUANPAN_FAIL;
    }

    converger_a = t_converger_a->get_copy();
    converger_b = t_converger_b->get_copy();

    converger_a->set_domain(t_domain);
    converger_b->set_domain(t_domain);

    return converger_a->initialize() + converger_b->initialize();
}

unique_ptr<Converger> LogicAND::get_copy() { return make_unique<LogicAND>(*this); }

bool LogicAND::is_converged(const unsigned counter) {
    auto result_a = std::async([&] { return converger_a->is_converged(counter); });
    auto result_b = std::async([&] { return converger_b->is_converged(counter); });

    const auto logic_result = result_a.get() && result_b.get();
    set_conv_flag(logic_result);
    converger_a->set_conv_flag(logic_result);
    converger_b->set_conv_flag(logic_result);

    return get_conv_flag();
}

unique_ptr<Converger> LogicOR::get_copy() { return make_unique<LogicOR>(*this); }

bool LogicOR::is_converged(const unsigned counter) {
    auto result_a = std::async([&] { return converger_a->is_converged(counter); });
    auto result_b = std::async([&] { return converger_b->is_converged(counter); });

    const auto logic_result = result_a.get() || result_b.get();
    set_conv_flag(logic_result);
    converger_a->set_conv_flag(logic_result);
    converger_b->set_conv_flag(logic_result);

    return get_conv_flag();
}

unique_ptr<Converger> LogicXOR::get_copy() { return make_unique<LogicXOR>(*this); }

bool LogicXOR::is_converged(const unsigned counter) {
    auto result_a = std::async([&] { return converger_a->is_converged(counter); });
    auto result_b = std::async([&] { return converger_b->is_converged(counter); });

    const auto logic_result = result_a.get() != result_b.get();
    set_conv_flag(logic_result);
    converger_a->set_conv_flag(logic_result);
    converger_b->set_conv_flag(logic_result);

    return get_conv_flag();
}
