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

#include "GroupParser.h"
#include <Domain/DomainBase.h>
#include <Domain/Group/CustomNodeGroup.h>
#include <Domain/Group/ElementGroup.h>
#include <Domain/Group/GroupGroup.h>
#include <Domain/Group/NodeGroup.h>
#include <Toolbox/utility.h>

using std::vector;

void new_nodegroup(unique_ptr<Group>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<uword> pool;
    get_input(command, pool);

    return_obj = make_unique<NodeGroup>(tag, pool);
}

void new_customnodegroup(unique_ptr<Group>& return_obj, istringstream& command) {
    unsigned tag, expression;
    if(!get_input(command, tag, expression)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    return_obj = make_unique<CustomNodeGroup>(tag, expression);
}

void new_elementgroup(unique_ptr<Group>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<uword> pool;
    get_input(command, pool);

    return_obj = make_unique<ElementGroup>(tag, pool);
}

void new_generate(unique_ptr<Group>& return_obj, istringstream& command) {
    string type;
    if(!get_input(command, type)) {
        suanpan_error("A valid type is required.\n");
        return;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    int start, interval, end;
    if(!get_input(command, start)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    if(!get_input(command, interval)) {
        interval = 1;
        end = start;
    }
    else if(!get_input(command, end)) {
        end = interval;
        interval = end > start ? 1 : -1;
    }

    if(0 == interval) interval = 1;

    if(start == end) interval = 1;
    else if(start < end && interval < 0 || start > end && interval > 0) interval = -interval;

    vector<uword> tag_pool;

    tag_pool.reserve(std::max(1, (end - start) / interval + 1));

    while(start <= end) {
        tag_pool.emplace_back(start);
        start += interval;
    }

    return_obj = make_unique<NodeGroup>(tag, tag_pool);

    if(is_equal(type, "nodegroup")) return_obj = make_unique<NodeGroup>(tag, tag_pool);
    else if(is_equal(type, "elementgroup")) return_obj = make_unique<ElementGroup>(tag, tag_pool);
}

void new_generatebyrule(unique_ptr<Group>& return_obj, istringstream& command) {
    if(string type; !get_input(command, type) || !is_equal(type, "nodegroup")) {
        suanpan_error("A valid type is required.\n");
        return;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned dof;
    if(!get_input(command, dof)) {
        suanpan_error("A valid dof identifier is required.\n");
        return;
    }

    vector<double> pool;
    get_input(command, pool);

    return_obj = make_unique<NodeGroup>(tag, dof, pool);
}

void new_generatebyplane(unique_ptr<Group>& return_obj, istringstream& command) {
    if(string type; !get_input(command, type) || !is_equal(type, "nodegroup")) {
        suanpan_error("A valid type is required.\n");
        return;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<double> pool;
    get_input(command, pool);

    if(pool.empty()) return;

    return_obj = make_unique<NodeGroup>(tag, pool);
}

void new_generatebypoint(unique_ptr<Group>& return_obj, istringstream& command) {
    if(string type; !get_input(command, type) || !is_equal(type, "nodegroup")) {
        suanpan_error("A valid type is required.\n");
        return;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<double> pool;
    get_input(command, pool);

    if(pool.size() % 2 != 0) return;

    const auto size = static_cast<long long>(pool.size()) / 2;

    return_obj = make_unique<NodeGroup>(tag, vector(pool.begin(), pool.begin() + size), vector(pool.end() - size, pool.end()));
}

void new_groupgroup(unique_ptr<Group>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<uword> pool;
    get_input(command, pool);

    return_obj = make_unique<GroupGroup>(tag, pool);
}

int create_new_group(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string group_id;
    if(!get_input(command, group_id)) {
        suanpan_error("A valid group type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Group> new_group = nullptr;

    if(is_equal(group_id, "nodegroup")) new_nodegroup(new_group, command);
    else if(is_equal(group_id, "customnodegroup")) new_customnodegroup(new_group, command);
    else if(is_equal(group_id, "elementgroup")) new_elementgroup(new_group, command);
    else if(is_equal(group_id, "groupgroup")) new_groupgroup(new_group, command);
    else if(is_equal(group_id, "generate")) new_generate(new_group, command);
    else if(is_equal(group_id, "generatebyrule")) new_generatebyrule(new_group, command);
    else if(is_equal(group_id, "generatebypoint")) new_generatebypoint(new_group, command);
    else if(is_equal(group_id, "generatebyplane")) new_generatebyplane(new_group, command);

    if(new_group == nullptr || !domain->insert(std::move(new_group)))
        suanpan_error("Fail to create new group via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
