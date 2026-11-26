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

// ReSharper disable StringLiteralTypo
// ReSharper disable IdentifierTypo
#include "ElementParser.h"

#include <Domain/ExternalModule.h>
#include <Element/Element>
#include <Toolbox/TokenRegistry.h>

namespace {
    //
    // elements
    //

    void new_allman(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Allman>(tag, std::move(node_tag), material_tag, thickness);
    }

    template<unsigned which> void new_b21(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(!get_optional_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        if(0 == which) return_obj = std::make_unique<B21>(tag, std::move(node_tag), section_id, int_pt, is_true(nonlinear));
        else return_obj = std::make_unique<B21E>(tag, which, std::move(node_tag), section_id, int_pt, is_true(nonlinear));
    }

    void new_b21h(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(!get_optional_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<B21H>(tag, std::move(node_tag), section_id, elastic_length, is_true(nonlinear));
    }

    template<bool if_os> void new_b31(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(!get_optional_input(command, nonlinear)) {
            suanpan_error("A valid nlgeom switch is required.\n");
            return;
        }

        if(if_os) return_obj = std::make_unique<B31OS>(tag, std::move(node_tag), section_id, orientation_id, int_pt, is_true(nonlinear));
        else return_obj = std::make_unique<B31>(tag, std::move(node_tag), section_id, orientation_id, int_pt, is_true(nonlinear));
    }

    template<unsigned which> void new_nmb21(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(!get_optional_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        if(0 == which) return_obj = std::make_unique<NMB21>(tag, std::move(node_tag), section_id, is_true(nonlinear));
        else return_obj = std::make_unique<NMB21E>(tag, which, std::move(node_tag), section_id, is_true(nonlinear));
    }

    void new_nmb31(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(!get_optional_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<NMB31>(tag, std::move(node_tag), section_id, orientation_id, is_true(nonlinear));
    }

    void new_c3d20(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string reduced_scheme = "true";
        if(command.eof())
            suanpan_debug("Standard Irons 14-point integration scheme assumed.\n");
        else if(!get_input(command, reduced_scheme))
            suanpan_error("A valid reduced integration switch is required.\n");

        std::string nonlinear = "false";
        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        return_obj = std::make_unique<C3D20>(tag, std::move(node_tag), material_tag, is_true(reduced_scheme), is_true(nonlinear));
    }

    void new_c3d4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        return_obj = std::make_unique<C3D4>(tag, std::move(node_tag), material_tag, is_true(nonlinear));
    }

    template<bool force_reduced> void new_c3d8(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string reduced_scheme = "I";
        if(force_reduced) reduced_scheme = "R";
        else if(command.eof()) suanpan_debug("Standard 4-point Gaussian integration scheme assumed.\n");
        else if(!get_input(command, reduced_scheme)) {
            suanpan_error("A valid integration switch is required.\n");
            return;
        }

        std::string nonlinear = "false";
        if(command.eof()) suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        auto penalty = 1.;
        if(command.eof()) suanpan_debug("Unit penalty multiplier assumed.\n");
        else if(!get_input(command, penalty)) {
            suanpan_error("A valid penalty multiplier is required.\n");
            return;
        }

        return_obj = std::make_unique<C3D8>(tag, std::move(node_tag), material_tag, penalty, suanpan::to_upper(reduced_scheme[0]), is_true(nonlinear));
    }

    void new_c3d8i(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<C3D8I>(tag, std::move(node_tag), material_tag);
    }

    void new_cax3(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        return_obj = std::make_unique<CAX3>(tag, std::move(node_tag), material_tag, is_true(nonlinear));
    }

    void new_cax4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<CAX4>(tag, std::move(node_tag), material_tag, false);
    }

    void new_cax8(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<CAX8>(tag, std::move(node_tag), material_tag, false);
    }

    void new_contact2d(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Contact2D>(tag, master_tag, slave_tag, alpha);
    }

    void new_contact3d(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Contact3D>(tag, master_tag, slave_tag, alpha);
    }

    void new_cp3(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        return_obj = std::make_unique<CP3>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
    }

    template<bool force_reduced> void new_cp4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string reduced_scheme = "false";
        if(force_reduced) reduced_scheme = "true";
        else if(command.eof()) suanpan_debug("Full integration assumed.\n");
        else if(!get_input(command, reduced_scheme)) {
            suanpan_error("A valid reduced scheme switch is required.\n");
            return;
        }

        std::string nonlinear = "false";
        if(command.eof()) suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        auto penalty = 1.;
        if(command.eof()) suanpan_debug("Unit penalty multiplier assumed.\n");
        else if(!get_input(command, penalty)) {
            suanpan_error("A valid penalty multiplier is required.\n");
            return;
        }

        return_obj = std::make_unique<CP4>(tag, std::move(node_tag), material_tag, thickness, penalty, is_true(reduced_scheme), is_true(nonlinear));
    }

    void new_cp4i(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<CP4I>(tag, std::move(node_tag), material_tag, thickness);
    }

    void new_cp5(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "N";
        if(!command.eof() && !get_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<CP5>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
    }

    void new_cp6(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(!get_optional_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<CP6>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
    }

    void new_cp7(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "N";
        if(!command.eof() && !get_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<CP7>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
    }

    void new_cp8(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string reduced_scheme = "N";
        if(!get_optional_input(command, reduced_scheme)) {
            suanpan_error("A valid reduced integration switch is required.\n");
            return;
        }

        std::string nonlinear = "N";
        if(!get_optional_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<CP8>(tag, std::move(node_tag), material_tag, thickness, is_true(reduced_scheme), is_true(nonlinear));
    }

    void new_cpe8(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string reduced_scheme = "N";
        if(!command.eof() && !get_input(command, reduced_scheme)) {
            suanpan_error("A valid reduced integration switch is required.\n");
            return;
        }

        std::string nonlinear = "N";
        if(!command.eof() && !get_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<CP8>(tag, std::move(node_tag), material_tag, 1., is_true(reduced_scheme), is_true(nonlinear));
    }

    void new_cpe8r(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "N";
        if(!command.eof() && !get_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<CP8>(tag, std::move(node_tag), material_tag, 1., true, is_true(nonlinear));
    }

    void new_cinp4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<CINP4>(tag, std::move(node_tag), material_tag, thickness);
    }

    void new_cin3d8(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<CIN3D8>(tag, std::move(node_tag), material_tag);
    }

    void new_csmt3(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<CSMT3>(tag, std::move(node_tag), material_tag, thickness, length);
    }

    void new_csmt6(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<CSMT6>(tag, std::move(node_tag), material_tag, thickness, length);
    }

    template<unsigned size> void new_csmq(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        if(4 == size) return_obj = std::make_unique<CSMQ4>(tag, std::move(node_tag), material_tag, thickness, length);
        else if(5 == size) return_obj = std::make_unique<CSMQ5>(tag, std::move(node_tag), material_tag, thickness, length);
        else if(6 == size) return_obj = std::make_unique<CSMQ6>(tag, std::move(node_tag), material_tag, thickness, length);
        else if(7 == size) return_obj = std::make_unique<CSMQ7>(tag, std::move(node_tag), material_tag, thickness, length);
        else if(8 == size) return_obj = std::make_unique<CSMQ8>(tag, std::move(node_tag), material_tag, thickness, length);
    }

    template<unsigned dimension> void new_damper01(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Damper01>(tag, std::move(node_tag), damper_tag, dimension);
    }

    template<unsigned dimension> void new_damper02(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string use_matrix = "true";
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

        return_obj = std::make_unique<Damper02>(tag, std::move(node_tag), damper_tag, spring_tag, is_true(use_matrix), proceed, beta, dimension);
    }

    template<unsigned dimension> void new_damper05(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Damper05>(tag, std::move(node_tag), damper_tag, dimension);
    }

    template<bool monolithic> void new_dc3d4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<DC3D4>(tag, std::move(node_tag), material_tag, length, rate, monolithic);
    }

    template<bool monolithic> void new_dc3d8(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<DC3D8>(tag, std::move(node_tag), material_tag, length, rate, monolithic);
    }

    void new_nonlocalc3d8(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<NonlocalC3D8>(tag, std::move(node_tag), material_tag);
    }

    template<bool monolithic> void new_dcp3(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<DCP3>(tag, std::move(node_tag), material_tag, length, rate, thickness, monolithic);
    }

    template<bool monolithic> void new_dcp4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<DCP4>(tag, std::move(node_tag), material_tag, length, rate, thickness, monolithic);
    }

    void new_nonlocalcp4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<NonlocalCP4>(tag, std::move(node_tag), material_tag, thickness);
    }

    void new_dkt3(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        auto num_ip = 3u;
        if(!get_optional_input(command, num_ip)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        return_obj = std::make_unique<DKT3>(tag, std::move(node_tag), material_tag, thickness, num_ip);
    }

    void new_dkt4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        auto num_ip = 3u;
        if(!get_optional_input(command, num_ip)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        return_obj = std::make_unique<DKT4>(tag, std::move(node_tag), material_tag, thickness, num_ip);
    }

    void new_dkts3(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        auto num_ip = 3u;
        if(!get_optional_input(command, num_ip)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        std::string nlgeom = "false";
        if(!get_optional_input(command, nlgeom)) {
            suanpan_error("A valid nlgeom switch is required.\n");
            return;
        }

        return_obj = std::make_unique<DKTS3>(tag, std::move(node_tag), material_tag, thickness, num_ip, is_true(nlgeom));
    }

    template<unsigned dof> void new_embedded(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Embedded<dof>>(tag, element_tag, node_tag, alpha);
    }

    void new_eb21(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        return_obj = std::make_unique<EB21>(tag, std::move(node_tag), area, moment_inertia, material_tag, is_true(nonlinear));
    }

    void new_eb31os(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        return_obj = std::make_unique<EB31OS>(tag, std::move(node_tag), std::move(property), orientation, is_true(nonlinear));
    }

    void new_f21(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        auto int_pt = 6u;
        if(!get_optional_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        std::string nonlinear = "false";
        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        return_obj = std::make_unique<F21>(tag, std::move(node_tag), section_id, int_pt, is_true(nonlinear));
    }

    void new_f21h(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        return_obj = std::make_unique<F21H>(tag, std::move(node_tag), section_id, elastic_length, is_true(nonlinear));
    }

    void new_f31(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "false";
        if(!command.eof() && !get_input(command, nonlinear)) {
            suanpan_error("A valid nonlinear geometry switch is required.\n");
            return;
        }

        return_obj = std::make_unique<F31>(tag, std::move(node_tag), section_id, orientation_id, int_pt, is_true(nonlinear));
    }

    void new_gcmq(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string int_scheme = "I";
        if(!command.eof() && !get_input(command, int_scheme))
            suanpan_error("A valid reduced scheme switch is required.\n");

        return_obj = std::make_unique<GCMQ>(tag, std::move(node_tag), material_tag, thickness, -1., suanpan::to_upper(int_scheme[0]));
    }

    template<bool simplified, char int_type> void new_gcmq(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        auto thickness = 1., objective_length = -1.;
        if(!get_optional_input(command, thickness)) {
            suanpan_error("A valid thickness is required.\n");
            return;
        }
        if(!get_optional_input(command, objective_length)) {
            suanpan_error("A valid length is required.\n");
            return;
        }

        if constexpr(simplified) return_obj = std::make_unique<SGCMQ>(tag, std::move(node_tag), material_tag, thickness, objective_length, int_type);
        else return_obj = std::make_unique<GCMQ>(tag, std::move(node_tag), material_tag, thickness, objective_length, int_type);
    }

    void new_sgcms(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nlgeom = "false";
        if(!get_optional_input(command, nlgeom)) {
            suanpan_error("A valid nlgeom switch is required.\n");
            return;
        }

        return_obj = std::make_unique<SGCMS>(tag, std::move(node_tag), material_tag, thickness, is_true(nlgeom));
    }

    void new_gq12(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<GQ12>(tag, std::move(node_tag), material_tag, thickness);
    }

    void new_joint(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
        uvec node_tag(2);
        if(!get_input(command, node_tag)) {
            suanpan_error("Two valid nodes are required.\n");
            return;
        }

        return_obj = std::make_unique<Joint>(tag, std::move(node_tag), get_remaining<uword>(command));
    }

    template<unsigned which> void new_mass(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        const auto dof_tag = get_remaining<uword>(command);

        if(2 == which && *std::ranges::max_element(dof_tag) > 3) {
            suanpan_error("At most three dofs are supported.\n");
            return;
        }
        if(3 == which && *std::ranges::max_element(dof_tag) > 6) {
            suanpan_error("At most six dofs are supported.\n");
            return;
        }

        if(2 == which) return_obj = std::make_unique<Mass2D>(tag, node, magnitude, uvec(dof_tag));
        else return_obj = std::make_unique<Mass3D>(tag, node, magnitude, uvec(dof_tag));
    }

    template<unsigned which> void new_masspoint(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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
            if(2 == which) return_obj = std::make_unique<MassPoint2D>(tag, node, translational_magnitude);
            else return_obj = std::make_unique<MassPoint3D>(tag, node, translational_magnitude);
        }
        else {
            double rotational_magnitude;
            if(!get_input(command, rotational_magnitude)) {
                suanpan_error("A valid rotational magnitude is required.\n");
                return;
            }
            if(2 == which) return_obj = std::make_unique<MassPoint2D>(tag, node, translational_magnitude, rotational_magnitude);
            else return_obj = std::make_unique<MassPoint3D>(tag, node, translational_magnitude, rotational_magnitude);
        }
    }

    void new_mindlin(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Mindlin>(tag, std::move(node_tag), material_tag, thickness, num_ip);
    }

    void new_mvlem(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        const auto [B, H, R, CT, ST] = get_remaining<double, double, double, uword, uword>(command);

        return_obj = std::make_unique<MVLEM>(tag, std::move(node_tag), B, H, R, CT, ST, shear_tag, c_height);
    }

    template<unsigned node> void new_pcpedc(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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
        if(!get_input(command, alpha, n, k)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        if(4 == node) return_obj = std::make_unique<PCPE4DC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n, k);
        else if(8 == node) return_obj = std::make_unique<PCPE8DC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n, k);
    }

    template<unsigned node> void new_pcpeuc(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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
        if(!get_input(command, alpha, n)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        if(4 == node) return_obj = std::make_unique<PCPE4UC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n);
        else if(8 == node) return_obj = std::make_unique<PCPE8UC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n);
    }

    void new_ps(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<PS>(tag, std::move(node_tag), material_tag, thickness);
    }

    void new_qe2(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<QE2>(tag, std::move(node_tag), material_tag, thickness);
    }

    void new_s4(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nlgeom = "false";
        if(!get_optional_input(command, nlgeom)) {
            suanpan_error("A valid nlgeom switch is required.\n");
            return;
        }

        return_obj = std::make_unique<S4>(tag, std::move(node_tag), material_tag, thickness, is_true(nlgeom));
    }

    template<typename T> void new_singlesection(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<T>(tag, node, section_tag);
    }

    template<unsigned DIM> void new_sphericalparticle(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
        unsigned node;
        if(!get_input(command, node)) {
            suanpan_error("A valid node tag is required.\n");
            return;
        }

        const auto pool = get_remaining<double>(command);
        if(pool.size() != 5 && pool.size() != 6) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        if constexpr(2u == DIM) {
            if(5 == pool.size())
                return_obj = std::make_unique<SphericalParticle2D>(tag, node, pool[0], pool[1], pool[2], pool[3], pool[4]);
            else
                return_obj = std::make_unique<InertialSphericalParticle2D>(tag, node, pool[0], pool[1], pool[2], pool[3], pool[4], pool[5]);
        }
        else {
            if(5 == pool.size())
                return_obj = std::make_unique<SphericalParticle3D>(tag, node, pool[0], pool[1], pool[2], pool[3], pool[4]);
            else
                return_obj = std::make_unique<InertialSphericalParticle3D>(tag, node, pool[0], pool[1], pool[2], pool[3], pool[4], pool[5]);
        }
    }
    void new_spring01(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Spring01>(tag, std::move(node_tag), material_tag);
    }

    void new_spring02(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        return_obj = std::make_unique<Spring02>(tag, std::move(node_tag), material_tag);
    }

    void new_t2d2(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "N", update_area = "N", log_strain = "N";
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

        return_obj = std::make_unique<T2D2>(tag, std::move(node_tag), material_tag, area, is_true(nonlinear), is_true(update_area), is_true(log_strain), rigidity);
    }

    void new_t2d2s(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "N", log_strain = "N";

        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        if(command.eof())
            suanpan_debug("Engineering strain assumed.\n");
        else if(!get_input(command, log_strain))
            suanpan_error("A valid switch to indicate if to use engineering strain is required.\n");

        return_obj = std::make_unique<T2D2S>(tag, std::move(node_tag), section_tag, is_true(nonlinear), is_true(log_strain));
    }

    void new_t3d2(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "N", update_area = "N", log_strain = "N";

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

        return_obj = std::make_unique<T3D2>(tag, std::move(node_tag), material_tag, area, is_true(nonlinear), is_true(update_area), is_true(log_strain));
    }

    void new_t3d2s(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        std::string nonlinear = "N", log_strain = "N";

        if(command.eof())
            suanpan_debug("Linear geometry assumed.\n");
        else if(!get_input(command, nonlinear))
            suanpan_error("A valid nonlinear geometry switch is required.\n");

        if(command.eof())
            suanpan_debug("Engineering strain assumed.\n");
        else if(!get_input(command, log_strain))
            suanpan_error("A valid switch to indicate if to use engineering strain is required.\n");

        return_obj = std::make_unique<T3D2S>(tag, std::move(node_tag), section_tag, is_true(nonlinear), is_true(log_strain));
    }

    void new_tie(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
        double magnitude, penalty;
        if(!get_input(command, magnitude) || !get_input(command, penalty)) {
            suanpan_error("A valid magnitude is required.\n");
            return;
        }

        const auto [node_tag, dof_tag, weight_tag] = get_remaining<uword, uword, double>(command);

        return_obj = std::make_unique<Tie>(tag, node_tag, dof_tag, weight_tag, magnitude, penalty);
    }

    template<unsigned S> void new_translationconnector(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        if(2u == S) return_obj = std::make_unique<TranslationConnector2D>(tag, std::move(nodes), penalty);
        else return_obj = std::make_unique<TranslationConnector3D>(tag, std::move(nodes), penalty);
    }

    void new_patchquad(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
        std::vector<uword> node_tag;
        std::vector<double> knot_x, knot_y;
        auto material_tag = -1;
        auto thickness = 1.;

        while(!command.eof()) {
            if(std::string parameter; get_input(command, parameter)) {
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

        return_obj = std::make_unique<PatchQuad>(tag, std::move(knot_x), std::move(knot_y), std::move(node_tag), static_cast<unsigned>(material_tag), thickness);
    }

    void new_patchcube(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
        std::vector<uword> node_tag;
        std::vector<double> knot_x, knot_y, knot_z;
        auto material_tag = -1;

        while(!command.eof()) {
            if(std::string parameter; get_input(command, parameter)) {
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

        return_obj = std::make_unique<PatchCube>(tag, std::move(knot_x), std::move(knot_y), std::move(knot_z), std::move(node_tag), static_cast<unsigned>(material_tag));
    }

    void new_global_mass(unique_ptr<Element>& return_obj, std::istringstream& command, const unsigned tag) {
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

        const auto dof_tag = get_remaining<uword>(command);

        if(dof_tag.empty()) {
            suanpan_error("At least one dof is required.\n");
            return;
        }

        if(*std::ranges::max_element(dof_tag) > 6) {
            suanpan_error("At most six dofs are supported.\n");
            return;
        }

        return_obj = std::make_unique<Mass3D>(tag, node, magnitude, uvec(dof_tag));
    }

    //
    // modifiers
    //

    void new_lumped_simple(unique_ptr<Modifier>& return_obj, std::istringstream& command, const unsigned tag) {
        return_obj = std::make_unique<LumpedSimple>(tag, get_remaining<uword>(command));
    }

    void new_lumped_scale(unique_ptr<Modifier>& return_obj, std::istringstream& command, const unsigned tag) {
        return_obj = std::make_unique<LumpedScale>(tag, get_remaining<uword>(command));
    }

    void new_rayleigh_modifier(unique_ptr<Modifier>& return_obj, std::istringstream& command, const unsigned tag) {
        double a, b, c, d;
        if(!get_input(command, a, b, c, d)) {
            suanpan_error("A valid value is required.\n");
            return;
        }

        return_obj = std::make_unique<Rayleigh>(tag, a, b, c, d, get_remaining<uword>(command));
    }

    void new_elemental_lee(unique_ptr<Modifier>& return_obj, std::istringstream& command, const unsigned tag) {
        double damping_ratio = .02;
        if(!command.eof() && !get_input(command, damping_ratio)) {
            suanpan_error("A valid damping ratio is required.\n");
            return;
        }

        return_obj = std::make_unique<ElementalLee>(tag, damping_ratio, get_remaining<uword>(command));
    }

    void new_linear_viscosity(unique_ptr<Modifier>& return_obj, std::istringstream& command, const unsigned tag) {
        double mu;
        if(!get_input(command, mu)) {
            suanpan_error("A valid value is required.\n");
            return;
        }

        return_obj = std::make_unique<LinearViscosity>(tag, mu, get_remaining<uword>(command));
    }

    void new_elemental_modal(unique_ptr<Modifier>& return_obj, std::istringstream& command, const unsigned tag) {
        double a, b;
        if(!get_input(command, a, b)) {
            suanpan_error("A valid value is required.\n");
            return;
        }

        return_obj = std::make_unique<ElementalModal>(tag, a, b, get_remaining<uword>(command));
    }

    template<bool group> void new_elemental_nonviscous(unique_ptr<Modifier>& return_obj, std::istringstream& command, const unsigned tag) {
        unsigned ele_tag;
        if(!get_input(command, ele_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        const auto [m_r, m_i, s_r, s_i] = get_remaining<double, double, double, double>(command);

        const vec m_imag{m_i}, s_imag{s_i};
        if(accu(m_imag) + accu(s_imag) > 1E-10) {
            suanpan_error("Parameters should be conjugate pairs.\n");
            return;
        }

        auto m = cx_vec{vec{m_r}, m_imag}, s = cx_vec{vec{s_r}, s_imag};

        if(const auto sum = accu(m % exp(-1E8 * s)); sum.real() * sum.real() + sum.imag() * sum.imag() > 1E-10) {
            suanpan_error("The provided kernel does not converge to zero.\n");
            return;
        }

        if constexpr(group) return_obj = std::make_unique<ElementalNonviscousGroup>(tag, std::move(m), std::move(s), ele_tag);
        else return_obj = std::make_unique<ElementalNonviscous>(tag, std::move(m), std::move(s), uvec{ele_tag});
    }

    //
    // orientations
    //

    template<typename T> void new_orientation_object(unique_ptr<Orientation>& return_obj, std::istringstream& command, const unsigned tag) {
        vec xyz(3);
        if(!get_input(command, xyz)) {
            suanpan_error("A valid vector is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag, std::move(xyz));
    }
} // namespace

int create_new_element(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string element_id;
    if(!get_input(command, element_id)) {
        suanpan_error("A valid element type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Element> new_element = nullptr;

    static const TokenRegistry<void (*)(unique_ptr<Element>&, std::istringstream&, unsigned)> registry{
        {"Allman", new_allman},
        {"B21", new_b21<0>},
        {"B21EH", new_b21<2>},
        {"B21EL", new_b21<1>},
        {"B21H", new_b21h},
        {"B31", new_b31<false>},
        {"B31OS", new_b31<true>},
        {"C3D4", new_c3d4},
        {"C3D8", new_c3d8<false>},
        {"C3D8I", new_c3d8i},
        {"C3D8R", new_c3d8<true>},
        {"C3D20", new_c3d20},
        {"CAX3", new_cax3},
        {"CAX4", new_cax4},
        {"CAX8", new_cax8},
        {"CIN3D8", new_cin3d8},
        {"CINP4", new_cinp4},
        {"Contact2D", new_contact2d},
        {"Contact3D", new_contact3d},
        {"CP3", new_cp3},
        {"CP4", new_cp4<false>},
        {"CP4I", new_cp4i},
        {"CP4R", new_cp4<true>},
        {"CP5", new_cp5},
        {"CP6", new_cp6},
        {"CP7", new_cp7},
        {"CP8", new_cp8},
        {"CPE8", new_cpe8},
        {"CPE8R", new_cpe8r},
        {"CPS8", new_cp8},
        {"CSMQ4", new_csmq<4>},
        {"CSMQ5", new_csmq<5>},
        {"CSMQ6", new_csmq<6>},
        {"CSMQ7", new_csmq<7>},
        {"CSMQ8", new_csmq<8>},
        {"CSMT3", new_csmt3},
        {"CSMT6", new_csmt6},
        {"Damper01", new_damper01<2>},
        {"Damper02", new_damper02<2>},
        {"Damper03", new_damper01<3>},
        {"Damper04", new_damper02<3>},
        {"Damper05", new_damper05<2>},
        {"Damper06", new_damper05<3>},
        {"DC3D4", new_dc3d4<false>},
        {"DC3D4M", new_dc3d4<true>},
        {"DC3D8", new_dc3d8<false>},
        {"DC3D8M", new_dc3d8<true>},
        {"DCP3", new_dcp3<false>},
        {"DCP3M", new_dcp3<true>},
        {"DCP4", new_dcp4<false>},
        {"DCP4M", new_dcp4<true>},
        {"DKT3", new_dkt3},
        {"DKT4", new_dkt4},
        {"DKTS3", new_dkts3},
        {"EB21", new_eb21},
        {"EB31OS", new_eb31os},
        {"Embedded2D", new_embedded<2u>},
        {"Embedded3D", new_embedded<3u>},
        {"F21", new_f21},
        {"F21H", new_f21h},
        {"F31", new_f31},
        {"GCMQ", new_gcmq},
        {"GCMQG", new_gcmq<false, 'G'>},
        {"GCMQI", new_gcmq<false, 'I'>},
        {"GCMQL", new_gcmq<false, 'L'>},
        {"GQ12", new_gq12},
        {"Joint", new_joint},
        {"Mass", new_mass<3>},
        {"Mass2D", new_mass<2>},
        {"Mass3D", new_mass<3>},
        {"MassPoint2D", new_masspoint<2>},
        {"MassPoint3D", new_masspoint<3>},
        {"Mindlin", new_mindlin},
        {"MVLEM", new_mvlem},
        {"NMB21", new_nmb21<0>},
        {"NMB21EH", new_nmb21<2>},
        {"NMB21EL", new_nmb21<1>},
        {"NMB31", new_nmb31},
        {"NonlocalC3D8", new_nonlocalc3d8},
        {"NonlocalCP4", new_nonlocalcp4},
        {"PatchCube", new_patchcube},
        {"PatchQuad", new_patchquad},
        {"PCPE4DC", new_pcpedc<4>},
        {"PCPE4UC", new_pcpeuc<4>},
        {"PCPE8DC", new_pcpedc<8>},
        {"PCPE8UC", new_pcpeuc<8>},
        {"PS", new_ps},
        {"QE2", new_qe2},
        {"S4", new_s4},
        {"SGCMQG", new_gcmq<true, 'G'>},
        {"SGCMQI", new_gcmq<true, 'I'>},
        {"SGCMQL", new_gcmq<true, 'L'>},
        {"SGCMS", new_sgcms},
        {"SingleSection2D", new_singlesection<SingleSection2D>},
        {"SingleSection3D", new_singlesection<SingleSection3D>},
        {"SphericalParticle2D", new_sphericalparticle<2u>},
        {"SphericalParticle3D", new_sphericalparticle<3u>},
        {"Spring01", new_spring01},
        {"Spring02", new_spring02},
        {"T2D2", new_t2d2},
        {"T2D2S", new_t2d2s},
        {"T3D2", new_t3d2},
        {"T3D2S", new_t3d2s},
        {"Tie", new_tie},
        {"TranslationConnector2D", new_translationconnector<2u>},
        {"TranslationConnector3D", new_translationconnector<3u>},
    };

    if(const auto factory = registry.fetch(element_id); nullptr != factory) factory(new_element, command, tag);
    else external_module::object(new_element, domain, element_id, command, tag);

    if(nullptr == new_element || !domain->insert(std::move(new_element)))
        suanpan_error("Fail to create new element via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_mass(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Element> new_mass = nullptr;

    new_global_mass(new_mass, command, tag);

    if(nullptr == new_mass || !domain->insert(std::move(new_mass)))
        suanpan_error("Fail to create new mass via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_modifier(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string modifier_type;
    if(!get_input(command, modifier_type)) {
        suanpan_error("A valid modifier type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Modifier> new_modifier = nullptr;

    static const TokenRegistry<void (*)(unique_ptr<Modifier>&, std::istringstream&, unsigned)> registry{
        {"LumpedSimple", new_lumped_simple},
        {"LumpedScale", new_lumped_scale},
        {"Rayleigh", new_rayleigh_modifier},
        {"ElementalLee", new_elemental_lee},
        {"LinearViscosity", new_linear_viscosity},
        {"ElementalModal", new_elemental_modal},
        {"ElementalNonviscous", new_elemental_nonviscous<false>},
        {"ElementalNonviscousGroup", new_elemental_nonviscous<true>},
    };

    if(const auto factory = registry.fetch(modifier_type); nullptr != factory) factory(new_modifier, command, tag);
    else external_module::object(new_modifier, domain, modifier_type, command, tag);

    if(nullptr == new_modifier || !domain->insert(std::move(new_modifier)))
        suanpan_error("Fail to create new modifier via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_orientation(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string file_type;
    if(!get_input(command, file_type)) {
        suanpan_error("A valid type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Orientation> new_orientation = nullptr;

    static const TokenRegistry<void (*)(unique_ptr<Orientation>&, std::istringstream&, unsigned)> registry{
        {"B3DL", new_orientation_object<B3DL>},
        {"B3DC", new_orientation_object<B3DC>},
        {"B3DOSL", new_orientation_object<B3DOSL>},
        {"B3DOSC", new_orientation_object<B3DOSC>},
    };

    if(const auto factory = registry.fetch(file_type); nullptr != factory) factory(new_orientation, command, tag);

    if(nullptr == new_orientation || !domain->insert(std::move(new_orientation)))
        suanpan_error("Fail to create new orientation via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
