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

#include "StepParser.h"
#include <Domain/DomainBase.h>
#include <Step/Step>
#include <Toolbox/utility.h>
#include <Solver/Integrator/Integrator.h>

int create_new_step(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string step_type;
    if(!get_input(command, step_type)) {
        SP_E("A valid step type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        SP_E("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(step_type, "Frequency")) {
        auto eigen_number = 1;
        if(!command.eof() && !get_input(command, eigen_number)) {
            SP_E("A valid number of eigenvalues is required.\n");
            return SUANPAN_SUCCESS;
        }

        char type = 's';
        if(!get_optional_input(command, type)) {
            SP_E("A valid eigenvalue type is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<Frequency>(tag, eigen_number, suanpan::to_upper(type)))) domain->set_current_step_tag(tag);
        else
            SP_E("Cannot create new step.\n");
    }
    else if(is_equal(step_type, "Buckling") || is_equal(step_type, "Buckle")) {
        if(domain->insert(make_shared<Buckle>(tag))) domain->set_current_step_tag(tag);
        else
            SP_E("Cannot create new step.\n");
    }
    else if(is_equal(step_type, "Optimization") || is_equal(step_type, "Optimisation")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            SP_E("A valid time period is required.\n");
            return SUANPAN_SUCCESS;
        }
        if(domain->insert(make_shared<Optimization>(tag, time))) domain->set_current_step_tag(tag);
        else
            SP_E("Cannot create new step.\n");
    }
    else if(is_equal(step_type, "Static")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            SP_E("A valid time period is required.\n");
            return SUANPAN_SUCCESS;
        }
        if(domain->insert(make_shared<Static>(tag, time))) domain->set_current_step_tag(tag);
        else
            SP_E("Cannot create new step.\n");
    }
    else if(is_equal(step_type, "Dynamic") || is_equal(step_type, "ImplicitDynamic")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            SP_E("A valid time period is required.\n");
            return SUANPAN_SUCCESS;
        }
        if(domain->insert(make_shared<Dynamic>(tag, time, IntegratorType::Implicit))) domain->set_current_step_tag(tag);
        else
            SP_E("Cannot create new step.\n");
    }
    else if(is_equal(step_type, "ExplicitDynamic")) {
        auto time = 1.;
        if(!command.eof() && !get_input(command, time)) {
            SP_E("A valid time period is required.\n");
            return SUANPAN_SUCCESS;
        }
        if(domain->insert(make_shared<Dynamic>(tag, time, IntegratorType::Explicit))) domain->set_current_step_tag(tag);
        else
            SP_E("Cannot create new step.\n");
    }
    else if(is_equal(step_type, "ArcLength")) {
        unsigned node;
        if(!get_input(command, node)) {
            SP_E("A valid node tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        unsigned dof;
        if(!get_input(command, dof)) {
            SP_E("A valid dof identifier is required.\n");
            return SUANPAN_SUCCESS;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            SP_E("A valid magnitude is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<ArcLength>(tag, node, dof, magnitude))) domain->set_current_step_tag(tag);
        else
            SP_E("Fail to create new step via \"{}\".\n", command.str());
    }
    else
        SP_E("Cannot identify step type.\n");

    return SUANPAN_SUCCESS;
}
