/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "LoadParser.h"
#include <Domain/DomainBase.h>
#include <Load/Load>
#include <Toolbox/utility.h>

int create_new_acceleration(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_acceleration() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_acceleration() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_acceleration() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_acceleration() needs dof.\n");
        return SUANPAN_SUCCESS;
    }

    uword node_id;
    vector<uword> node_pool;
    while(get_input(command, node_id)) node_pool.emplace_back(node_id);

    if(const auto step_tag = domain->get_current_step_tag(); !domain->insert(make_shared<NodalAcceleration>(load_id, step_tag, magnitude, uvec(node_pool), dof_id, amplitude_id))) suanpan_error("create_new_acceleration() fails to create new load.\n");

    return SUANPAN_SUCCESS;
}

int create_new_amplitude(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string amplitude_type;
    if(!get_input(command, amplitude_type)) {
        suanpan_error("create_new_amplitude() needs a valid amplitude type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_amplitude() needs a valid amplitude type.\n");
        return SUANPAN_SUCCESS;
    }

    if(const auto step_tag = domain->get_current_step_tag(); is_equal(amplitude_type, "Constant")) domain->insert(make_shared<Constant>(tag, step_tag));
    else if(is_equal(amplitude_type, "Ramp")) domain->insert(make_shared<Ramp>(tag, step_tag));
    else if(is_equal(amplitude_type, "Tabular")) {
        string file_name;
        if(!get_input(command, file_name)) {
            suanpan_error("create_new_amplitude() needs a valid file.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(make_shared<Tabular>(tag, std::move(file_name), step_tag));
    }
    else if(is_equal(amplitude_type, "Decay")) {
        double A;
        if(!get_input(command, A)) {
            suanpan_error("create_new_amplitude() needs a A.\n");
            return SUANPAN_SUCCESS;
        }
        double TD;
        if(!get_input(command, TD)) {
            suanpan_error("create_new_amplitude() needs a TD.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(make_shared<Decay>(tag, A, TD, step_tag));
    }
    else if(is_equal(amplitude_type, "Linear")) {
        double A;
        if(!get_input(command, A)) {
            suanpan_error("create_new_amplitude() needs a slope.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(make_shared<Linear>(tag, A, step_tag));
    }
    else if(is_equal(amplitude_type, "Combine")) {
        vector<uword> tag_pool;
        uword t_tag;
        while(get_input(command, t_tag)) tag_pool.emplace_back(t_tag);
        domain->insert(make_shared<Combine>(tag, uvec(tag_pool), step_tag));
    }
    else if(is_equal(amplitude_type, "Modulated") || is_equal(amplitude_type, "Sine") || is_equal(amplitude_type, "Cosine")) {
        double W;
        if(!get_input(command, W)) {
            suanpan_error("create_new_amplitude() needs a period/amplitude.\n");
            return SUANPAN_SUCCESS;
        }

        double amp;
        vector<double> A;
        while(get_input(command, amp)) A.emplace_back(amp);

        if(is_equal(amplitude_type, "Modulated")) domain->insert(make_shared<Modulated>(tag, W, std::move(A), step_tag));
        else if(is_equal(amplitude_type, "Sine")) domain->insert(make_shared<Sine>(tag, W, std::move(A), step_tag));
        else if(is_equal(amplitude_type, "Cosine")) domain->insert(make_shared<Cosine>(tag, W, std::move(A), step_tag));
    }
    else if(is_equal(amplitude_type, "NZStrongMotion")) {
        string name;
        if(!get_input(command, name)) {
            suanpan_error("create_new_amplitude() needs a name.\n");
            return SUANPAN_SUCCESS;
        }

        domain->insert(make_shared<NZStrongMotion>(tag, name.c_str(), step_tag));
    }

    return SUANPAN_SUCCESS;
}

int create_new_bodyforce(const shared_ptr<DomainBase>& domain, istringstream& command, const bool flag) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_bodyforce() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_bodyforce() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_bodyforce() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_bodyforce() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned element;
    vector<uword> element_tag;
    while(get_input(command, element)) element_tag.push_back(element);

    if(flag) { if(!domain->insert(make_shared<GroupBodyForce>(load_id, domain->get_current_step_tag(), magnitude, uvec(element_tag), dof_id, amplitude_id))) suanpan_error("create_new_bodyforce() fails to create new load.\n"); }
    else if(!domain->insert(make_shared<BodyForce>(load_id, domain->get_current_step_tag(), magnitude, uvec(element_tag), dof_id, amplitude_id))) suanpan_error("create_new_bodyforce() fails to create new load.\n");

    return SUANPAN_SUCCESS;
}

int create_new_cload(const shared_ptr<DomainBase>& domain, istringstream& command, const bool flag) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_cload() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_cload() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_cload() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_cload() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(flag) { if(!domain->insert(make_shared<GroupNodalForce>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_cload() fails to create new load.\n"); }
    else { if(!domain->insert(make_shared<NodalForce>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_cload() fails to create new load.\n"); }

    return SUANPAN_SUCCESS;
}

int create_new_lineudl(const shared_ptr<DomainBase>& domain, istringstream& command, const unsigned dimension) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_lineudl() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_lineudl() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_lineudl() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_lineudl() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(2 == dimension) { if(!domain->insert(make_shared<LineUDL2D>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_lineudl() fails to create new load.\n"); }
    else { if(!domain->insert(make_shared<LineUDL3D>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_lineudl() fails to create new load.\n"); }

    return SUANPAN_SUCCESS;
}

int create_new_displacement(const shared_ptr<DomainBase>& domain, istringstream& command, const bool flag) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_displacement() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_displacement() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_displacement() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_displacement() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(const auto step_tag = domain->get_current_step_tag(); flag) { if(!domain->insert(make_shared<GroupNodalDisplacement>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_displacement() fails to create new load.\n"); }
    else { if(!domain->insert(make_shared<NodalDisplacement>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_displacement() fails to create new load.\n"); }

    return SUANPAN_SUCCESS;
}

int create_new_supportmotion(const shared_ptr<DomainBase>& domain, istringstream& command, const unsigned flag) {
    unsigned load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("create_new_supportmotion() needs a tag.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned amplitude_id;
    if(!get_input(command, amplitude_id)) {
        suanpan_error("create_new_supportmotion() needs a valid amplitude tag.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("create_new_supportmotion() needs load magnitude.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("create_new_supportmotion() needs a valid DoF.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    vector<uword> node_tag;
    while(get_input(command, node)) node_tag.push_back(node);

    if(const auto step_tag = domain->get_current_step_tag(); 0 == flag) { if(!domain->insert(make_shared<SupportDisplacement>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_supportmotion() fails to create new load.\n"); }
    else if(1 == flag) { if(!domain->insert(make_shared<SupportVelocity>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_supportmotion() fails to create new load.\n"); }
    else { if(!domain->insert(make_shared<SupportAcceleration>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_supportmotion() fails to create new load.\n"); }

    return SUANPAN_SUCCESS;
}
