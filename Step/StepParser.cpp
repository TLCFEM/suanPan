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

#include "StepParser.h"

#include <Domain/DomainBase.h>
#include <Step/Step>
#include <Toolbox/TokenRegistry.h>

namespace {
    void new_frequency(unique_ptr<Step>& return_obj, std::istringstream& command, const unsigned tag) {
        auto eigen_number = 1;
        if(!command.eof() && !get_input(command, eigen_number)) {
            suanpan_error("A valid number of eigenvalues is required.\n");
            return;
        }

        char type = 's';
        if(!get_optional_input(command, type)) {
            suanpan_error("A valid eigenvalue type is required.\n");
            return;
        }

        return_obj = std::make_unique<Frequency>(tag, eigen_number, suanpan::to_upper(type));
    }

    void new_buckle(unique_ptr<Step>& return_obj, std::istringstream& command, const unsigned tag) {
        return_obj = std::make_unique<Buckle>(tag);
    }

    template<typename T> void new_time_step(unique_ptr<Step>& return_obj, std::istringstream& command, const unsigned tag) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            suanpan_error("A valid time period is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag, time);
    }

    template<Integrator::Type type> void new_dynamic(unique_ptr<Step>& return_obj, std::istringstream& command, const unsigned tag) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            suanpan_error("A valid time period is required.\n");
            return;
        }

        return_obj = std::make_unique<Dynamic>(tag, time, type);
    }

    void new_arclength(unique_ptr<Step>& return_obj, std::istringstream& command, const unsigned tag) {
        return_obj = std::make_unique<ArcLength>(tag);
    }
} // namespace

int create_new_step(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string step_type;
    if(!get_input(command, step_type)) {
        suanpan_error("A valid step type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Step> new_step = nullptr;

    static const TokenRegistry<void (*)(unique_ptr<Step>&, std::istringstream&, unsigned)> registry{
        {"Frequency", new_frequency},
        {"Buckling", new_buckle},
        {"Buckle", new_buckle},
        {"Optimization", new_time_step<Optimization>},
        {"Optimisation", new_time_step<Optimization>},
        {"Static", new_time_step<Static>},
        {"Dynamic", new_dynamic<Integrator::Type::Implicit>},
        {"ImplicitDynamic", new_dynamic<Integrator::Type::Implicit>},
        {"ExplicitDynamic", new_dynamic<Integrator::Type::Explicit>},
        {"ArcLength", new_arclength},
    };

    if(const auto factory = registry.fetch(step_type); nullptr != factory) factory(new_step, command, tag);

    if(nullptr == new_step || !domain->insert(std::move(new_step))) {
        suanpan_error("Fail to create new step via \"{}\".\n", command.str());
        return SUANPAN_SUCCESS;
    }

    domain->set_current_step_tag(tag);

    return SUANPAN_SUCCESS;
}
