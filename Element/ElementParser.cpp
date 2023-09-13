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

#include "ElementParser.h"
#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Element/Element>
#include <Toolbox/utility.h>

using std::vector;

void new_allman(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(3);
    if(!get_input(command, node_tag)) {
        suanpan_error("Three valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_debug("Unit thickness assumed.\n");
    else if(!get_input(command, thickness))
        suanpan_error("A valid thickness is required.\n");

    return_obj = make_unique<Allman>(tag, std::move(node_tag), material_tag, thickness);
}

void new_b21(unique_ptr<Element>& return_obj, istringstream& command, const unsigned which) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!get_optional_input(command, int_pt)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    if(0 == which) return_obj = make_unique<B21>(tag, std::move(node_tag), section_id, int_pt, is_true(nonlinear));
    else return_obj = make_unique<B21E>(tag, which, std::move(node_tag), section_id, int_pt, is_true(nonlinear));
}

void new_b21h(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    auto elastic_length = .2;
    if(!command.eof() && !get_input(command, elastic_length)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<B21H>(tag, std::move(node_tag), section_id, elastic_length, is_true(nonlinear));
}

void new_b31(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    unsigned orientation_id;
    if(!get_input(command, orientation_id)) {
        suanpan_error("A valid orientation tag is required.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!get_optional_input(command, int_pt)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("A valid nlgeom switch is required.\n");
        return;
    }

    return_obj = make_unique<B31>(tag, std::move(node_tag), section_id, orientation_id, int_pt, is_true(nonlinear));
}

void new_nmb21(unique_ptr<Element>& return_obj, istringstream& command, const unsigned which) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    if(0 == which) return_obj = make_unique<NMB21>(tag, std::move(node_tag), section_id, is_true(nonlinear));
    else return_obj = make_unique<NMB21E>(tag, which, std::move(node_tag), section_id, is_true(nonlinear));
}

void new_nmb31(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    unsigned orientation_id;
    if(!get_input(command, orientation_id)) {
        suanpan_error("A valid orientation tag is required.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<NMB31>(tag, std::move(node_tag), section_id, orientation_id, is_true(nonlinear));
}

void new_c3d20(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(20);
    if(!get_input(command, node_tag)) {
        suanpan_error("Twenty valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    string reduced_scheme = "true";
    if(command.eof())
        suanpan_debug("Standard Irons 14-point integration scheme assumed.\n");
    else if(!get_input(command, reduced_scheme))
        suanpan_error("A valid reduced integration switch is required.\n");

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<C3D20>(tag, std::move(node_tag), material_tag, is_true(reduced_scheme), is_true(nonlinear));
}

void new_c3d4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<C3D4>(tag, std::move(node_tag), material_tag, is_true(nonlinear));
}

void new_c3d8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Eight valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    string reduced_scheme = "I";
    if(command.eof())
        suanpan_debug("Standard 4-point integration scheme assumed.\n");
    else if(!get_input(command, reduced_scheme))
        suanpan_error("A valid reduced integration switch is required.\n");

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<C3D8>(tag, std::move(node_tag), material_tag, suanpan::to_upper(reduced_scheme[0]), is_true(nonlinear));
}

void new_c3d8r(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Eight valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<C3D8>(tag, std::move(node_tag), material_tag, 'R', is_true(nonlinear));
}

void new_c3d8i(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Eight valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    return_obj = make_unique<C3D8I>(tag, std::move(node_tag), material_tag);
}

void new_cax3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(3);
    if(!get_input(command, node_tag)) {
        suanpan_error("Three valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<CAX3>(tag, std::move(node_tag), material_tag, is_true(nonlinear));
}

void new_cax4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    return_obj = make_unique<CAX4>(tag, std::move(node_tag), material_tag, false);
}

void new_cax8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Eight valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    return_obj = make_unique<CAX8>(tag, std::move(node_tag), material_tag, false);
}

void new_contact2d(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned master_tag, slave_tag;
    if(!get_input(command, master_tag)) {
        suanpan_error("A valid master group tag is required.\n");
        return;
    }
    if(!get_input(command, slave_tag)) {
        suanpan_error("A valid slave group tag is required.\n");
        return;
    }

    auto alpha = 1E6;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("A valid multiplier is required.\n");
        return;
    }

    return_obj = make_unique<Contact2D>(tag, master_tag, slave_tag, alpha);
}

void new_contact3d(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned master_tag, slave_tag;
    if(!get_input(command, master_tag)) {
        suanpan_error("A valid master group tag is required.\n");
        return;
    }
    if(!get_input(command, slave_tag)) {
        suanpan_error("A valid slave group tag is required.\n");
        return;
    }

    auto alpha = 1E6;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("A valid multiplier is required.\n");
        return;
    }

    return_obj = make_unique<Contact3D>(tag, master_tag, slave_tag, alpha);
}

void new_cp3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(3);
    if(!get_input(command, node_tag)) {
        suanpan_error("Three valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<CP3>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
}

void new_cp4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string reduced_scheme = "N";
    if(command.eof())
        suanpan_debug("Full integration assumed.\n");
    else if(!get_input(command, reduced_scheme)) {
        suanpan_error("A valid reduced scheme switch is required.\n");
        return;
    }

    string nonlinear = "N";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<CP4>(tag, std::move(node_tag), material_tag, thickness, is_true(reduced_scheme), is_true(nonlinear));
}

void new_cp4i(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    return_obj = make_unique<CP4I>(tag, std::move(node_tag), material_tag, thickness);
}

void new_cp4r(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<CP4>(tag, std::move(node_tag), material_tag, thickness, true, is_true(nonlinear));
}

void new_cp5(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(5);
    if(!get_input(command, node_tag)) {
        suanpan_error("Five valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<CP5>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
}

void new_cp6(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(6);
    if(!get_input(command, node_tag)) {
        suanpan_error("Six valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<CP6>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
}

void new_cp7(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(7);
    if(!get_input(command, node_tag)) {
        suanpan_error("Five valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<CP7>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
}

void new_cp8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Eight valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string reduced_scheme = "N";
    if(!get_optional_input(command, reduced_scheme)) {
        suanpan_error("A valid reduced integration switch is required.\n");
        return;
    }

    string nonlinear = "N";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<CP8>(tag, std::move(node_tag), material_tag, thickness, is_true(reduced_scheme), is_true(nonlinear));
}

void new_cpe8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Eight valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    string reduced_scheme = "N";
    if(!command.eof() && !get_input(command, reduced_scheme)) {
        suanpan_error("A valid reduced integration switch is required.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<CP8>(tag, std::move(node_tag), material_tag, 1., is_true(reduced_scheme), is_true(nonlinear));
}

void new_cpe8r(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Eight valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<CP8>(tag, std::move(node_tag), material_tag, 1., true, is_true(nonlinear));
}

void new_cinp4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    return_obj = make_unique<CINP4>(tag, std::move(node_tag), material_tag, thickness);
}

void new_cin3d8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Eight valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    return_obj = make_unique<CIN3D8>(tag, std::move(node_tag), material_tag);
}

void new_csmt3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(3);
    if(!get_input(command, node_tag)) {
        suanpan_error("Three valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    auto length = -1.;
    if(!get_optional_input(command, length)) {
        suanpan_error("A valid length is required.\n");
        return;
    }

    return_obj = make_unique<CSMT3>(tag, std::move(node_tag), material_tag, thickness, length);
}

void new_csmt6(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(6);
    if(!get_input(command, node_tag)) {
        suanpan_error("Six valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    auto length = -1.;
    if(!get_optional_input(command, length)) {
        suanpan_error("A valid length is required.\n");
        return;
    }

    return_obj = make_unique<CSMT6>(tag, std::move(node_tag), material_tag, thickness, length);
}

void new_csmq(unique_ptr<Element>& return_obj, istringstream& command, const unsigned size) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(size);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    auto length = -1.;
    if(!get_optional_input(command, length)) {
        suanpan_error("A valid length is required.\n");
        return;
    }

    if(4 == size) return_obj = make_unique<CSMQ4>(tag, std::move(node_tag), material_tag, thickness, length);
    else if(5 == size) return_obj = make_unique<CSMQ5>(tag, std::move(node_tag), material_tag, thickness, length);
    else if(6 == size) return_obj = make_unique<CSMQ6>(tag, std::move(node_tag), material_tag, thickness, length);
    else if(7 == size) return_obj = make_unique<CSMQ7>(tag, std::move(node_tag), material_tag, thickness, length);
    else if(8 == size) return_obj = make_unique<CSMQ8>(tag, std::move(node_tag), material_tag, thickness, length);
}

void new_damper01(unique_ptr<Element>& return_obj, istringstream& command, const unsigned dimension) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned damper_tag;
    if(!get_input(command, damper_tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    return_obj = make_unique<Damper01>(tag, std::move(node_tag), damper_tag, dimension);
}

void new_damper02(unique_ptr<Element>& return_obj, istringstream& command, const unsigned dimension) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned damper_tag;
    if(!get_input(command, damper_tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned spring_tag;
    if(!get_input(command, spring_tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    string use_matrix = "true";
    if(!command.eof() && !get_input(command, use_matrix)) {
        suanpan_error("A valid switch is required.\n");
        return;
    }

    unsigned proceed = 0;
    if(!command.eof() && !get_input(command, proceed)) {
        suanpan_error("A valid proceed switch is required.\n");
        return;
    }

    auto beta = .5;
    if(!command.eof() && !get_input(command, beta)) {
        suanpan_error("A valid beta value is required.\n");
        return;
    }

    return_obj = make_unique<Damper02>(tag, std::move(node_tag), damper_tag, spring_tag, is_true(use_matrix), proceed, beta, dimension);
}

void new_damper05(unique_ptr<Element>& return_obj, istringstream& command, const unsigned dimension) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned damper_tag;
    if(!get_input(command, damper_tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    return_obj = make_unique<Damper05>(tag, std::move(node_tag), damper_tag, dimension);
}

void new_dc3d4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double length;
    if(!get_input(command, length)) {
        suanpan_error("A valid characteristic length is required.\n");
        return;
    }

    double rate;
    if(!get_input(command, rate)) {
        suanpan_error("A valid energy release rate is required.\n");
        return;
    }

    return_obj = make_unique<DC3D4>(tag, std::move(node_tag), material_tag, length, rate);
}

void new_dc3d8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(8);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double length;
    if(!get_input(command, length)) {
        suanpan_error("A valid characteristic length is required.\n");
        return;
    }

    double rate;
    if(!get_input(command, rate)) {
        suanpan_error("A valid energy release rate is required.\n");
        return;
    }

    return_obj = make_unique<DC3D8>(tag, std::move(node_tag), material_tag, length, rate);
}

void new_dcp3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(3);
    if(!get_input(command, node_tag)) {
        suanpan_error("Three valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double length;
    if(!get_input(command, length)) {
        suanpan_error("A valid characteristic length is required.\n");
        return;
    }

    double rate;
    if(!get_input(command, rate)) {
        suanpan_error("A valid energy release rate is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    return_obj = make_unique<DCP3>(tag, std::move(node_tag), material_tag, length, rate, thickness);
}

void new_dcp4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double length;
    if(!get_input(command, length)) {
        suanpan_error("A valid characteristic length is required.\n");
        return;
    }

    double rate;
    if(!get_input(command, rate)) {
        suanpan_error("A valid energy release rate is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    return_obj = make_unique<DCP4>(tag, std::move(node_tag), material_tag, length, rate, thickness);
}

void new_dkt3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(3);
    if(!get_input(command, node_tag)) {
        suanpan_error("Three valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double thickness;
    if(!get_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    unsigned num_ip = 3;
    if(!get_optional_input(command, num_ip)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    return_obj = make_unique<DKT3>(tag, std::move(node_tag), material_tag, thickness, num_ip);
}

void new_dkt4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double thickness;
    if(!get_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    unsigned num_ip = 3;
    if(!get_optional_input(command, num_ip)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    return_obj = make_unique<DKT4>(tag, std::move(node_tag), material_tag, thickness, num_ip);
}

void new_dkts3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(3);
    if(!get_input(command, node_tag)) {
        suanpan_error("Three valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double thickness;
    if(!get_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    unsigned num_ip = 3;
    if(!get_optional_input(command, num_ip)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    string nlgeom = "false";
    if(!get_optional_input(command, nlgeom)) {
        suanpan_error("A valid nlgeom switch is required.\n");
        return;
    }

    return_obj = make_unique<DKTS3>(tag, std::move(node_tag), material_tag, thickness, num_ip, is_true(nlgeom));
}

void new_embedded(unique_ptr<Element>& return_obj, istringstream& command, const unsigned dof) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned element_tag;
    if(!get_input(command, element_tag)) {
        suanpan_error("A valid element tag is required.\n");
        return;
    }

    unsigned node_tag;
    if(!get_input(command, node_tag)) {
        suanpan_error("A valid node tag is required.\n");
        return;
    }

    auto alpha = 1E4;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("A valid alpha is required.\n");
        return;
    }

    if(2 == dof) return_obj = make_unique<Embedded2D>(tag, element_tag, node_tag, alpha);
    else return_obj = make_unique<Embedded3D>(tag, element_tag, node_tag, alpha);
}

void new_eb21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    double area;
    if(!get_input(command, area)) {
        suanpan_error("A valid area is required.\n");
        return;
    }

    double moment_inertia;
    if(!get_input(command, moment_inertia)) {
        suanpan_error("A valid moment of inertia is required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<EB21>(tag, std::move(node_tag), area, moment_inertia, material_tag, is_true(nonlinear));
}

void new_eb31os(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    vec property(7);
    if(!get_input(command, property)) {
        suanpan_error("A valid section/material property is required.\n");
        return;
    }

    unsigned orientation;
    if(!get_input(command, orientation)) {
        suanpan_error("A valid orientation tag is required.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<EB31OS>(tag, std::move(node_tag), std::move(property), orientation, is_true(nonlinear));
}

void new_f21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!get_optional_input(command, int_pt)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<F21>(tag, std::move(node_tag), section_id, int_pt, !!nonlinear);
}

void new_f21h(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    auto elastic_length = .2;
    if(!command.eof() && !get_input(command, elastic_length)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    return_obj = make_unique<F21H>(tag, std::move(node_tag), section_id, elastic_length, !!nonlinear);
}

void new_f31(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    unsigned orientation_id;
    if(!get_input(command, orientation_id)) {
        suanpan_error("A valid orientation tag is required.\n");
        return;
    }

    unsigned int_pt = 5;
    if(!command.eof() && !get_input(command, int_pt)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    string nonlinear = "false";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("A valid nonlinear geometry switch is required.\n");
        return;
    }

    return_obj = make_unique<F31>(tag, std::move(node_tag), section_id, orientation_id, int_pt, is_true(nonlinear));
}

void new_gcmq(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_debug("Unit thickness assumed.\n");
    else if(!get_input(command, thickness))
        suanpan_error("A valid thickness is required.\n");

    string int_scheme = "I";
    if(!command.eof() && !get_input(command, int_scheme))
        suanpan_error("A valid reduced scheme switch is required.\n");

    return_obj = make_unique<GCMQ>(tag, std::move(node_tag), material_tag, thickness, suanpan::to_upper(int_scheme[0]));
}

void new_gcmq(unique_ptr<Element>& return_obj, istringstream& command, const char int_type) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_debug("Unit thickness assumed.\n");
    else if(!get_input(command, thickness))
        suanpan_error("A valid thickness is required.\n");

    return_obj = make_unique<GCMQ>(tag, std::move(node_tag), material_tag, thickness, int_type);
}

void new_sgcmq(unique_ptr<Element>& return_obj, istringstream& command, const char int_type) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_debug("Unit thickness assumed.\n");
    else if(!get_input(command, thickness))
        suanpan_error("A valid thickness is required.\n");

    return_obj = make_unique<SGCMQ>(tag, std::move(node_tag), material_tag, thickness, int_type);
}

void new_sgcms(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && !get_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string nlgeom = "false";
    if(!get_optional_input(command, nlgeom)) {
        suanpan_error("A valid nlgeom switch is required.\n");
        return;
    }

    return_obj = make_unique<SGCMS>(tag, std::move(node_tag), material_tag, thickness, is_true(nlgeom));
}

void new_gq12(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_debug("Unit thickness assumed.\n");
    else if(!get_input(command, thickness))
        suanpan_error("A valid thickness is required.\n");

    return_obj = make_unique<GQ12>(tag, std::move(node_tag), material_tag, thickness);
}

void new_joint(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    vector<uword> pool;
    uword material_tag;
    while(!command.eof() && get_input(command, material_tag)) pool.emplace_back(material_tag);

    return_obj = make_unique<Joint>(tag, std::move(node_tag), std::move(pool));
}

void new_mass(unique_ptr<Element>& return_obj, istringstream& command, const unsigned which) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("A valid node tag is required.\n");
        return;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("A valid magnitude is required.\n");
        return;
    }

    unsigned dof;
    vector<uword> dof_tag;
    while(!command.eof() && get_input(command, dof)) dof_tag.push_back(dof);

    if(2 == which && *std::max_element(dof_tag.cbegin(), dof_tag.cend()) > 3) {
        suanpan_error("At most three dofs are supported.\n");
        return;
    }
    if(3 == which && *std::max_element(dof_tag.cbegin(), dof_tag.cend()) > 6) {
        suanpan_error("At most six dofs are supported.\n");
        return;
    }

    if(2 == which) return_obj = make_unique<Mass2D>(tag, node, magnitude, uvec(dof_tag));
    else return_obj = make_unique<Mass3D>(tag, node, magnitude, uvec(dof_tag));
}

void new_masspoint(unique_ptr<Element>& return_obj, istringstream& command, const unsigned which) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("A valid node tag is required.\n");
        return;
    }

    double translational_magnitude;
    if(!get_input(command, translational_magnitude)) {
        suanpan_error("A valid translational magnitude is required.\n");
        return;
    }

    if(command.eof()) {
        if(2 == which) return_obj = make_unique<MassPoint2D>(tag, node, translational_magnitude);
        else return_obj = make_unique<MassPoint3D>(tag, node, translational_magnitude);
    }
    else {
        double rotational_magnitude;
        if(!get_input(command, rotational_magnitude)) {
            suanpan_error("A valid rotational magnitude is required.\n");
            return;
        }
        if(2 == which) return_obj = make_unique<MassPoint2D>(tag, node, translational_magnitude, rotational_magnitude);
        else return_obj = make_unique<MassPoint3D>(tag, node, translational_magnitude, rotational_magnitude);
    }
}

void new_mindlin(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double thickness;
    if(!get_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    unsigned num_ip = 5;
    if(!command.eof() && !get_input(command, num_ip)) {
        suanpan_error("A valid number of integration points is required.\n");
        return;
    }

    return_obj = make_unique<Mindlin>(tag, std::move(node_tag), material_tag, thickness, num_ip);
}

void new_mvlem(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned shear_tag;
    if(!get_input(command, shear_tag)) {
        suanpan_error("A valid shear material is required.\n");
        return;
    }

    double c_height;
    if(!get_input(command, c_height)) {
        suanpan_error("A valid c is required.\n");
        return;
    }

    vector<double> B, H, R;
    vector<uword> CT, ST;
    while(!command.eof()) {
        double t_value;
        uword t_tag;
        if(!get_input(command, t_value)) {
            suanpan_error("A valid fibre width is required.\n");
            return;
        }
        B.emplace_back(t_value);
        if(!get_input(command, t_value)) {
            suanpan_error("A valid fibre thickness is required.\n");
            return;
        }
        H.emplace_back(t_value);
        if(!get_input(command, t_value)) {
            suanpan_error("A valid fibre reinforcement ratio is required.\n");
            return;
        }
        R.emplace_back(t_value);
        if(!get_input(command, t_tag)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }
        CT.emplace_back(t_tag);
        if(!get_input(command, t_tag)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }
        ST.emplace_back(t_tag);
    }

    return_obj = make_unique<MVLEM>(tag, std::move(node_tag), B, H, R, uvec(CT), uvec(ST), shear_tag, c_height);
}

void new_pcpedc(unique_ptr<Element>& return_obj, istringstream& command, const unsigned node) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(node);
    if(!get_input(command, node_tag)) {
        suanpan_error("{} valid nodes are required.\n", node);
        return;
    }

    unsigned solid_tag;
    if(!get_input(command, solid_tag)) {
        suanpan_error("A valid material tag for solid phase is required.\n");
        return;
    }
    unsigned fluid_tag;
    if(!get_input(command, fluid_tag)) {
        suanpan_error("A valid material tag for fluid phase is required.\n");
        return;
    }

    double alpha, n, k;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("A valid alpha is required.\n");
        return;
    }
    if(!get_optional_input(command, n)) {
        suanpan_error("A valid porosity is required.\n");
        return;
    }
    if(!get_optional_input(command, k)) {
        suanpan_error("A valid permeability is required.\n");
        return;
    }

    if(4 == node) return_obj = make_unique<PCPE4DC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n, k);
    else if(8 == node) return_obj = make_unique<PCPE8DC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n, k);
}

void new_pcpeuc(unique_ptr<Element>& return_obj, istringstream& command, const unsigned node) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(node);
    if(!get_input(command, node_tag)) {
        suanpan_error("{} valid nodes are required.\n", node);
        return;
    }

    unsigned solid_tag;
    if(!get_input(command, solid_tag)) {
        suanpan_error("A valid material tag for solid phase is required.\n");
        return;
    }
    unsigned fluid_tag;
    if(!get_input(command, fluid_tag)) {
        suanpan_error("A valid material tag for fluid phase is required.\n");
        return;
    }

    double alpha, n;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("A valid alpha is required.\n");
        return;
    }
    if(!get_optional_input(command, n)) {
        suanpan_error("A valid porosity is required.\n");
        return;
    }

    if(4 == node) return_obj = make_unique<PCPE4UC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n);
    else if(8 == node) return_obj = make_unique<PCPE8UC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n);
}

void new_ps(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_debug("Unit thickness assumed.\n");
    else if(!get_input(command, thickness))
        suanpan_error("A valid thickness is required.\n");

    return_obj = make_unique<PS>(tag, std::move(node_tag), material_tag, thickness);
}

void new_qe2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof())
        suanpan_debug("Unit thickness assumed.\n");
    else if(!get_input(command, thickness))
        suanpan_error("A valid thickness is required.\n");

    return_obj = make_unique<QE2>(tag, std::move(node_tag), material_tag, thickness);
}

void new_s4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(4);
    if(!get_input(command, node_tag)) {
        suanpan_error("Four valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && !get_input(command, thickness)) {
        suanpan_error("A valid thickness is required.\n");
        return;
    }

    string nlgeom = "false";
    if(!get_optional_input(command, nlgeom)) {
        suanpan_error("A valid nlgeom switch is required.\n");
        return;
    }

    return_obj = make_unique<S4>(tag, std::move(node_tag), material_tag, thickness, is_true(nlgeom));
}

void new_singlesection2d(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("A valid node tag is required.\n");
        return;
    }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    return_obj = make_unique<SingleSection2D>(tag, node, section_tag);
}

void new_singlesection3d(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("A valid node tag is required.\n");
        return;
    }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("A valid section tag is required.\n");
        return;
    }

    return_obj = make_unique<SingleSection3D>(tag, node, section_tag);
}

void new_spring01(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    return_obj = make_unique<Spring01>(tag, std::move(node_tag), material_tag);
}

void new_spring02(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    return_obj = make_unique<Spring02>(tag, std::move(node_tag), material_tag);
}

void new_t2d2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double area;
    if(!get_input(command, area)) {
        suanpan_error("A valid area is required.\n");
        return;
    }

    string nonlinear = "N", update_area = "N", log_strain = "N";
    double rigidity = -1.;

    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    if(command.eof())
        suanpan_debug("Constant area assumed.\n");
    else if(!get_input(command, update_area))
        suanpan_error("A valid area switch is required.\n");

    if(command.eof())
        suanpan_debug("Engineering strain assumed.\n");
    else if(!get_input(command, log_strain))
        suanpan_error("A valid engineering strain switch is required.\n");

    if(command.eof())
        suanpan_debug("No buckling check.\n");
    else if(!get_input(command, rigidity))
        suanpan_error("A valid flexural rigidity is required.\n");

    return_obj = make_unique<T2D2>(tag, std::move(node_tag), material_tag, area, is_true(nonlinear), is_true(update_area), is_true(log_strain), rigidity);
}

void new_t2d2s(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("A valid material/section tag is required.\n");
        return;
    }

    string nonlinear = "N", log_strain = "N";

    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    if(command.eof())
        suanpan_debug("Engineering strain assumed.\n");
    else if(!get_input(command, log_strain))
        suanpan_error("A valid switch to indicate if to use engineering strain is required.\n");

    return_obj = make_unique<T2D2S>(tag, std::move(node_tag), section_tag, is_true(nonlinear), is_true(log_strain));
}

void new_t3d2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    double area;
    if(!get_input(command, area)) {
        suanpan_error("A valid area is required.\n");
        return;
    }

    string nonlinear = "N", update_area = "N", log_strain = "N";

    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    if(command.eof())
        suanpan_debug("Constant area assumed.\n");
    else if(!get_input(command, update_area))
        suanpan_error("A valid area switch is required.\n");

    if(command.eof())
        suanpan_debug("Engineering strain assumed.\n");
    else if(!get_input(command, log_strain))
        suanpan_error("A valid engineering strain switch is required.\n");

    return_obj = make_unique<T3D2>(tag, std::move(node_tag), material_tag, area, is_true(nonlinear), is_true(update_area), is_true(log_strain));
}

void new_t3d2s(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec node_tag(2);
    if(!get_input(command, node_tag)) {
        suanpan_error("Two valid nodes are required.\n");
        return;
    }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("A valid material/section tag is required.\n");
        return;
    }

    string nonlinear = "N", log_strain = "N";

    if(command.eof())
        suanpan_debug("Linear geometry assumed.\n");
    else if(!get_input(command, nonlinear))
        suanpan_error("A valid nonlinear geometry switch is required.\n");

    if(command.eof())
        suanpan_debug("Engineering strain assumed.\n");
    else if(!get_input(command, log_strain))
        suanpan_error("A valid switch to indicate if to use engineering strain is required.\n");

    return_obj = make_unique<T3D2S>(tag, std::move(node_tag), section_tag, is_true(nonlinear), is_true(log_strain));
}

void new_tie(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    double magnitude, penalty;
    if(!get_input(command, magnitude) || !get_input(command, penalty)) {
        suanpan_error("A valid magnitude is required.\n");
        return;
    }

    vector<uword> node_tag, dof_tag;
    vector<double> weight_tag;
    while(!command.eof()) {
        double weight;
        uword dof;
        uword node;
        if(!get_input(command, node) || !get_input(command, dof) || !get_input(command, weight)) return;
        node_tag.emplace_back(node);
        dof_tag.emplace_back(dof);
        weight_tag.emplace_back(weight);
    }

    return_obj = make_unique<Tie>(tag, node_tag, dof_tag, weight_tag, magnitude, penalty);
}

void new_translationconnector(unique_ptr<Element>& return_obj, istringstream& command, const unsigned S) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    uvec nodes(3);
    if(!get_input(command, nodes)) {
        suanpan_error("Three valid nodes are required.\n");
        return;
    }

    double penalty;
    if(!get_input(command, penalty)) {
        suanpan_error("A valid penalty is required.\n");
        return;
    }

    if(2u == S) return_obj = make_unique<TranslationConnector2D>(tag, std::move(nodes), penalty);
    else return_obj = make_unique<TranslationConnector3D>(tag, std::move(nodes), penalty);
}

void new_patchquad(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<uword> node_tag;
    vector<double> knot_x, knot_y;
    auto material_tag = -1;
    auto thickness = 1.;

    while(!command.eof()) {
        if(string parameter; get_input(command, parameter)) {
            ignore_whitespace(command);
            if(uword node; is_equal(parameter, "-node"))
                while(command.peek() != '-' && get_input(command, node)) {
                    node_tag.emplace_back(node);
                    ignore_whitespace(command);
                }
            else if(double knot; is_equal(parameter, "-knotx"))
                while(command.peek() != '-' && get_input(command, knot)) {
                    knot_x.emplace_back(knot);
                    ignore_whitespace(command);
                }
            else if(is_equal(parameter, "-knoty"))
                while(command.peek() != '-' && get_input(command, knot)) {
                    knot_y.emplace_back(knot);
                    ignore_whitespace(command);
                }
            else if(is_equal(parameter, "-material")) {
                if(!get_input(command, material_tag)) {
                    suanpan_error("A valid material tag is required.\n");
                    return;
                }
                ignore_whitespace(command);
            }
            else if(is_equal(parameter, "-thickness")) {
                if(!get_input(command, thickness)) {
                    suanpan_error("A valid thickness is required.\n");
                    return;
                }
                ignore_whitespace(command);
            }
        }
    }

    if(node_tag.empty()) {
        suanpan_error("A valid node tags is required.\n");
        return;
    }
    if(knot_x.empty() || knot_y.empty()) {
        suanpan_error("A valid knot vector is required.\n");
        return;
    }
    if(-1 == material_tag) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    return_obj = make_unique<PatchQuad>(tag, std::move(knot_x), std::move(knot_y), std::move(node_tag), static_cast<unsigned>(material_tag), thickness);
}

void new_patchcube(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<uword> node_tag;
    vector<double> knot_x, knot_y, knot_z;
    auto material_tag = -1;

    while(!command.eof()) {
        if(string parameter; get_input(command, parameter)) {
            ignore_whitespace(command);
            if(uword node; is_equal(parameter, "-node"))
                while(command.peek() != '-' && get_input(command, node)) {
                    node_tag.emplace_back(node);
                    ignore_whitespace(command);
                }
            else if(double knot; is_equal(parameter, "-knotx"))
                while(command.peek() != '-' && get_input(command, knot)) {
                    knot_x.emplace_back(knot);
                    ignore_whitespace(command);
                }
            else if(is_equal(parameter, "-knoty"))
                while(command.peek() != '-' && get_input(command, knot)) {
                    knot_y.emplace_back(knot);
                    ignore_whitespace(command);
                }
            else if(is_equal(parameter, "-knotz"))
                while(command.peek() != '-' && get_input(command, knot)) {
                    knot_z.emplace_back(knot);
                    ignore_whitespace(command);
                }
            else if(is_equal(parameter, "-material")) {
                if(!get_input(command, material_tag)) {
                    suanpan_error("A valid material tag is required.\n");
                    return;
                }
                ignore_whitespace(command);
            }
        }
    }

    if(node_tag.empty()) {
        suanpan_error("A valid node tags is required.\n");
        return;
    }
    if(knot_x.empty() || knot_y.empty() || knot_z.empty()) {
        suanpan_error("A valid knot vector is required.\n");
        return;
    }
    if(-1 == material_tag) {
        suanpan_error("A valid material tag is required.\n");
        return;
    }

    return_obj = make_unique<PatchCube>(tag, std::move(knot_x), std::move(knot_y), std::move(knot_z), std::move(node_tag), static_cast<unsigned>(material_tag));
}

int create_new_mass(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("A valid node tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("A valid magnitude is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof;
    vector<uword> dof_tag;
    while(get_input(command, dof)) dof_tag.push_back(dof);

    if(*std::max_element(dof_tag.cbegin(), dof_tag.cend()) > 6) {
        suanpan_error("At most six dofs are supported.\n");
        return SUANPAN_SUCCESS;
    }

    domain->insert(make_shared<Mass3D>(tag, node, magnitude, uvec(dof_tag)));

    return SUANPAN_SUCCESS;
}

int create_new_modifier(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string modifier_type;
    if(!get_input(command, modifier_type)) {
        suanpan_error("A valid modifier type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Modifier> new_modifier = nullptr;

    auto get_element_pool = [&] {
        vector<uword> element_tag;
        unsigned e_tag;
        while(!command.eof() && get_input(command, e_tag)) element_tag.emplace_back(e_tag);

        return uvec(element_tag);
    };

    if(is_equal(modifier_type, "LumpedSimple")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        new_modifier = make_unique<LumpedSimple>(tag, get_element_pool());
    }
    else if(is_equal(modifier_type, "LumpedScale")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        new_modifier = make_unique<LumpedScale>(tag, get_element_pool());
    }
    else if(is_equal(modifier_type, "Rayleigh")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        double a, b, c, d;
        if(!get_input(command, a, b, c, d)) {
            suanpan_error("A valid value is required.\n");
            return SUANPAN_SUCCESS;
        }

        new_modifier = make_unique<Rayleigh>(tag, a, b, c, d, get_element_pool());
    }
    else if(is_equal(modifier_type, "LeeElementalDamping")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        double damping_ratio = -1.;
        if(!command.eof() && !get_input(command, damping_ratio)) {
            suanpan_error("A valid damping ratio is required.\n");
            return SUANPAN_SUCCESS;
        }

        new_modifier = make_unique<LeeElementalDamping>(tag, damping_ratio, get_element_pool());
    }
    else if(is_equal(modifier_type, "LinearViscosity")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        double mu;
        if(!get_input(command, mu)) {
            suanpan_error("A valid value is required.\n");
            return SUANPAN_SUCCESS;
        }

        new_modifier = make_unique<LinearViscosity>(tag, mu, get_element_pool());
    }
    else if(is_equal(modifier_type, "ElementalModal")) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        double a, b;
        if(!get_input(command, a, b)) {
            suanpan_error("A valid value is required.\n");
            return SUANPAN_SUCCESS;
        }

        new_modifier = make_unique<ElementalModal>(tag, a, b, get_element_pool());
    }
    else {
        // check if the library is already loaded
        auto code = false;
        for(const auto& I : domain->get_external_module_pool())
            if(is_equal(I->library_name, modifier_type) || I->locate_cpp_module(modifier_type)) {
                code = true;
                break;
            }

        // not loaded then try load it
        if(!code && domain->insert(make_shared<ExternalModule>(modifier_type))) code = true;

        // if loaded find corresponding function
        if(code)
            for(const auto& I : domain->get_external_module_pool()) {
                if(I->locate_cpp_module(modifier_type)) I->new_object(new_modifier, command);
                if(new_modifier != nullptr) break;
            }
    }

    if(nullptr == new_modifier || !domain->insert(std::move(new_modifier)))
        suanpan_error("Fail to create new modifier via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_orientation(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string file_type;
    if(!get_input(command, file_type)) {
        suanpan_error("A valid type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    vec xyz(3);
    if(!get_input(command, xyz)) {
        suanpan_error("A valid vector is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(file_type, "B3DL")) domain->insert(make_shared<B3DL>(tag, std::move(xyz)));
    else if(is_equal(file_type, "B3DC")) domain->insert(make_shared<B3DC>(tag, std::move(xyz)));
    else if(is_equal(file_type, "B3DOSL")) domain->insert(make_shared<B3DOSL>(tag, std::move(xyz)));
    else if(is_equal(file_type, "B3DOSC")) domain->insert(make_shared<B3DOSC>(tag, std::move(xyz)));

    return SUANPAN_SUCCESS;
}

int create_new_element(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string element_id;
    if(!get_input(command, element_id)) {
        suanpan_error("A valid element type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Element> new_element = nullptr;

    if(is_equal(element_id, "Allman")) new_allman(new_element, command);
    else if(is_equal(element_id, "B21")) new_b21(new_element, command, 0);
    else if(is_equal(element_id, "B21EL")) new_b21(new_element, command, 1);
    else if(is_equal(element_id, "B21EH")) new_b21(new_element, command, 2);
    else if(is_equal(element_id, "B21H")) new_b21h(new_element, command);
    else if(is_equal(element_id, "B31")) new_b31(new_element, command);
    else if(is_equal(element_id, "NMB21")) new_nmb21(new_element, command, 0);
    else if(is_equal(element_id, "NMB21EL")) new_nmb21(new_element, command, 1);
    else if(is_equal(element_id, "NMB21EH")) new_nmb21(new_element, command, 2);
    else if(is_equal(element_id, "NMB31")) new_nmb31(new_element, command);
    else if(is_equal(element_id, "C3D20")) new_c3d20(new_element, command);
    else if(is_equal(element_id, "C3D4")) new_c3d4(new_element, command);
    else if(is_equal(element_id, "C3D8")) new_c3d8(new_element, command);
    else if(is_equal(element_id, "C3D8R")) new_c3d8r(new_element, command);
    else if(is_equal(element_id, "C3D8I")) new_c3d8i(new_element, command);
    else if(is_equal(element_id, "CAX3")) new_cax3(new_element, command);
    else if(is_equal(element_id, "CAX4")) new_cax4(new_element, command);
    else if(is_equal(element_id, "CAX8")) new_cax8(new_element, command);
    else if(is_equal(element_id, "CIN3D8")) new_cin3d8(new_element, command);
    else if(is_equal(element_id, "CINP4")) new_cinp4(new_element, command);
    else if(is_equal(element_id, "Contact2D")) new_contact2d(new_element, command);
    else if(is_equal(element_id, "Contact3D")) new_contact3d(new_element, command);
    else if(is_equal(element_id, "CP3")) new_cp3(new_element, command);
    else if(is_equal(element_id, "CP4")) new_cp4(new_element, command);
    else if(is_equal(element_id, "PCPE4DC")) new_pcpedc(new_element, command, 4);
    else if(is_equal(element_id, "PCPE8DC")) new_pcpedc(new_element, command, 8);
    else if(is_equal(element_id, "PCPE4UC")) new_pcpeuc(new_element, command, 4);
    else if(is_equal(element_id, "PCPE8UC")) new_pcpeuc(new_element, command, 8);
    else if(is_equal(element_id, "CP4I")) new_cp4i(new_element, command);
    else if(is_equal(element_id, "CP4R")) new_cp4r(new_element, command);
    else if(is_equal(element_id, "CP5")) new_cp5(new_element, command);
    else if(is_equal(element_id, "CP6")) new_cp6(new_element, command);
    else if(is_equal(element_id, "CP7")) new_cp7(new_element, command);
    else if(is_equal(element_id, "CP8")) new_cp8(new_element, command);
    else if(is_equal(element_id, "CPE8")) new_cpe8(new_element, command);
    else if(is_equal(element_id, "CPE8R")) new_cpe8r(new_element, command);
    else if(is_equal(element_id, "CPS8")) new_cp8(new_element, command);
    else if(is_equal(element_id, "CSMT3")) new_csmt3(new_element, command);
    else if(is_equal(element_id, "CSMT6")) new_csmt6(new_element, command);
    else if(is_equal(element_id, "CSMQ4")) new_csmq(new_element, command, 4);
    else if(is_equal(element_id, "CSMQ5")) new_csmq(new_element, command, 5);
    else if(is_equal(element_id, "CSMQ6")) new_csmq(new_element, command, 6);
    else if(is_equal(element_id, "CSMQ7")) new_csmq(new_element, command, 7);
    else if(is_equal(element_id, "CSMQ8")) new_csmq(new_element, command, 8);
    else if(is_equal(element_id, "Damper01")) new_damper01(new_element, command, 2);
    else if(is_equal(element_id, "Damper02")) new_damper02(new_element, command, 2);
    else if(is_equal(element_id, "Damper03")) new_damper01(new_element, command, 3);
    else if(is_equal(element_id, "Damper04")) new_damper02(new_element, command, 3);
    else if(is_equal(element_id, "Damper05")) new_damper05(new_element, command, 2);
    else if(is_equal(element_id, "Damper06")) new_damper05(new_element, command, 3);
    else if(is_equal(element_id, "DC3D4")) new_dc3d4(new_element, command);
    else if(is_equal(element_id, "DC3D8")) new_dc3d8(new_element, command);
    else if(is_equal(element_id, "DCP3")) new_dcp3(new_element, command);
    else if(is_equal(element_id, "DCP4")) new_dcp4(new_element, command);
    else if(is_equal(element_id, "DKT3")) new_dkt3(new_element, command);
    else if(is_equal(element_id, "DKT4")) new_dkt4(new_element, command);
    else if(is_equal(element_id, "DKTS3")) new_dkts3(new_element, command);
    else if(is_equal(element_id, "Embedded2D")) new_embedded(new_element, command, 2);
    else if(is_equal(element_id, "Embedded3D")) new_embedded(new_element, command, 3);
    else if(is_equal(element_id, "EB21")) new_eb21(new_element, command);
    else if(is_equal(element_id, "EB31OS")) new_eb31os(new_element, command);
    else if(is_equal(element_id, "F21")) new_f21(new_element, command);
    else if(is_equal(element_id, "F21H")) new_f21h(new_element, command);
    else if(is_equal(element_id, "F31")) new_f31(new_element, command);
    else if(is_equal(element_id, "GCMQ")) new_gcmq(new_element, command);
    else if(is_equal(element_id, "GCMQG")) new_gcmq(new_element, command, 'G');
    else if(is_equal(element_id, "GCMQI")) new_gcmq(new_element, command, 'I');
    else if(is_equal(element_id, "GCMQL")) new_gcmq(new_element, command, 'L');
    else if(is_equal(element_id, "GQ12")) new_gq12(new_element, command);
    else if(is_equal(element_id, "Joint")) new_joint(new_element, command);
    else if(is_equal(element_id, "Mass2D")) new_mass(new_element, command, 2);
    else if(is_equal(element_id, "Mass3D")) new_mass(new_element, command, 3);
    else if(is_equal(element_id, "Mass")) new_mass(new_element, command, 3);
    else if(is_equal(element_id, "MassPoint2D")) new_masspoint(new_element, command, 2);
    else if(is_equal(element_id, "MassPoint3D")) new_masspoint(new_element, command, 3);
    else if(is_equal(element_id, "Mindlin")) new_mindlin(new_element, command);
    else if(is_equal(element_id, "MVLEM")) new_mvlem(new_element, command);
    else if(is_equal(element_id, "PS")) new_ps(new_element, command);
    else if(is_equal(element_id, "QE2")) new_qe2(new_element, command);
    else if(is_equal(element_id, "S4")) new_s4(new_element, command);
    else if(is_equal(element_id, "SGCMQG")) new_sgcmq(new_element, command, 'G');
    else if(is_equal(element_id, "SGCMQI")) new_sgcmq(new_element, command, 'I');
    else if(is_equal(element_id, "SGCMQL")) new_sgcmq(new_element, command, 'L');
    else if(is_equal(element_id, "SGCMS")) new_sgcms(new_element, command);
    else if(is_equal(element_id, "SingleSection2D")) new_singlesection2d(new_element, command);
    else if(is_equal(element_id, "SingleSection3D")) new_singlesection3d(new_element, command);
    else if(is_equal(element_id, "Spring01")) new_spring01(new_element, command);
    else if(is_equal(element_id, "Spring02")) new_spring02(new_element, command);
    else if(is_equal(element_id, "T2D2")) new_t2d2(new_element, command);
    else if(is_equal(element_id, "T2D2S")) new_t2d2s(new_element, command);
    else if(is_equal(element_id, "T3D2")) new_t3d2(new_element, command);
    else if(is_equal(element_id, "T3D2S")) new_t3d2s(new_element, command);
    else if(is_equal(element_id, "Tie")) new_tie(new_element, command);
    else if(is_equal(element_id, "TranslationConnector2D")) new_translationconnector(new_element, command, 2u);
    else if(is_equal(element_id, "TranslationConnector3D")) new_translationconnector(new_element, command, 3u);
    else if(is_equal(element_id, "PatchQuad")) new_patchquad(new_element, command);
    else if(is_equal(element_id, "PatchCube")) new_patchcube(new_element, command);
    else load::object(new_element, domain, element_id, command);

    if(new_element == nullptr || !domain->insert(std::move(new_element)))
        suanpan_error("Fail to create new element via \"{}\".\n", command.str());

    return 0;
}
