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

#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Element/Element>
#include <Toolbox/utility.h>

using std::vector;

int create_new_element(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string element_id;
    if(!get_input(command, element_id)) {
        suanpan_error("create_new_element() needs element type.\n");
        return 0;
    }

    unique_ptr<Element> new_element = nullptr;

    if(is_equal(element_id, "Allman")) new_allman(new_element, command);
    else if(is_equal(element_id, "B21")) new_b21(new_element, command);
    else if(is_equal(element_id, "B21EL")) new_b21e(new_element, command, 1);
    else if(is_equal(element_id, "B21EH")) new_b21e(new_element, command, 2);
    else if(is_equal(element_id, "B21H")) new_b21h(new_element, command);
    else if(is_equal(element_id, "B31")) new_b31(new_element, command);
    else if(is_equal(element_id, "NMB21")) new_nmb21(new_element, command);
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
    else if(is_equal(element_id, "Damper01")) new_damper01(new_element, command);
    else if(is_equal(element_id, "Damper02")) new_damper02(new_element, command);
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
    else if(is_equal(element_id, "F21")) new_f21(new_element, command);
    else if(is_equal(element_id, "F21H")) new_f21h(new_element, command);
    else if(is_equal(element_id, "F31")) new_f31(new_element, command);
    else if(is_equal(element_id, "GCMQ")) new_gcmq(new_element, command);
    else if(is_equal(element_id, "GCMQG")) new_gcmqg(new_element, command);
    else if(is_equal(element_id, "GCMQI")) new_gcmqi(new_element, command);
    else if(is_equal(element_id, "GCMQL")) new_gcmql(new_element, command);
    else if(is_equal(element_id, "GQ12")) new_gq12(new_element, command);
    else if(is_equal(element_id, "Joint")) new_joint(new_element, command);
    else if(is_equal(element_id, "Mass")) new_mass(new_element, command);
    else if(is_equal(element_id, "Mindlin")) new_mindlin(new_element, command);
    else if(is_equal(element_id, "MVLEM")) new_mvlem(new_element, command);
    else if(is_equal(element_id, "PS")) new_ps(new_element, command);
    else if(is_equal(element_id, "QE2")) new_qe2(new_element, command);
    else if(is_equal(element_id, "S4")) new_s4(new_element, command);
    else if(is_equal(element_id, "SGCMQG")) new_sgcmqg(new_element, command);
    else if(is_equal(element_id, "SGCMQI")) new_sgcmqi(new_element, command);
    else if(is_equal(element_id, "SGCMQL")) new_sgcmql(new_element, command);
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
    else if(is_equal(element_id, "PatchQuad")) new_patchquad(new_element, command);
    else if(is_equal(element_id, "PatchCube")) new_patchcube(new_element, command);
    else load::object(new_element, domain, element_id, command);

    if(new_element == nullptr || !domain->insert(std::move(new_element))) suanpan_error("create_new_element() fails to create new element.\n");

    return 0;
}

void new_allman(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_allman() needs a tag.\n");
        return;
    }

    uvec node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_allman() needs three valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_allman() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_allman() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_allman() needs a valid thickness.\n");

    return_obj = make_unique<Allman>(tag, std::move(node_tag), material_tag, thickness);
}

void new_b21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_b21() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_b21() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_b21() needs a valid section tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!get_optional_input(command, int_pt)) {
        suanpan_error("new_b21() needs a valid number of integration points.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("new_b21() needs a valid nonlinear geometry switch (0,1).\n");
        return;
    }

    return_obj = make_unique<B21>(tag, std::move(node_tag), section_id, int_pt, is_true(nonlinear));
}

void new_b21e(unique_ptr<Element>& return_obj, istringstream& command, const unsigned which) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_b21e() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_b21e() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_b21e() needs a valid section tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!get_optional_input(command, int_pt)) {
        suanpan_error("new_b21e() needs a valid number of integration points.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("new_b21e() needs a valid nonlinear geometry switch (0,1).\n");
        return;
    }

    return_obj = make_unique<B21E>(tag, which, std::move(node_tag), section_id, int_pt, is_true(nonlinear));
}

void new_b21h(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_b21h() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_b21h() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_b21h() needs a valid section tag.\n");
        return;
    }

    auto elastic_length = .2;
    if(!command.eof() && !get_input(command, elastic_length)) {
        suanpan_error("new_b21h() needs a valid number of integration points.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("new_b21h() needs a valid nonlinear geometry switch (0,1).\n");
        return;
    }

    return_obj = make_unique<B21H>(tag, std::move(node_tag), section_id, elastic_length, is_true(nonlinear));
}

void new_b31(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_b31() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_b31() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_b31() needs a valid section tag.\n");
        return;
    }

    unsigned orientation_id;
    if(!get_input(command, orientation_id)) {
        suanpan_error("new_b31() needs a valid orientation tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!get_optional_input(command, int_pt)) {
        suanpan_error("new_b31() needs a valid number of integration points.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("new_b31() needs a valid nlgeom switch.\n");
        return;
    }

    return_obj = make_unique<B31>(tag, std::move(node_tag), section_id, orientation_id, int_pt, is_true(nonlinear));
}

void new_nmb21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_nmb21() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_nmb21() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_nmb21() needs a valid section tag.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("new_nmb21() needs a valid nonlinear geometry switch.\n");
        return;
    }

    return_obj = make_unique<NMB21>(tag, std::move(node_tag), section_id, is_true(nonlinear));
}

void new_nmb31(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_nmb31() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_nmb31() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_nmb31() needs a valid section tag.\n");
        return;
    }

    unsigned orientation_id;
    if(!get_input(command, orientation_id)) {
        suanpan_error("new_nmb31() needs a valid orientation tag.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("new_nmb31() needs a valid nonlinear geometry switch.\n");
        return;
    }

    return_obj = make_unique<NMB31>(tag, std::move(node_tag), section_id, orientation_id, is_true(nonlinear));
}

void new_c3d20(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_c3d20() needs a valid tag.\n");
        return;
    }

    uvec node_tag(20);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_c3d20() needs twenty valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_c3d20() needs a valid material tag.\n");
        return;
    }

    string reduced_scheme = "true";
    if(command.eof()) suanpan_debug("new_c3d20() assumes standard Irons 14-point integration scheme.\n");
    else if(!get_input(command, reduced_scheme)) suanpan_error("new_c3d20() needs a valid reduced integration switch (0,1).\n");

    string nonlinear = "false";
    if(command.eof()) suanpan_debug("new_c3d20() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_c3d20() needs a valid nonlinear geometry switch (0,1).\n");

    return_obj = make_unique<C3D20>(tag, std::move(node_tag), material_tag, is_true(reduced_scheme), is_true(nonlinear));
}

void new_c3d4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_c3d4() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_c3d4() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_c3d4() needs a valid material tag.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof()) suanpan_debug("new_c3d4() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_c3d4() needs a valid nonlinear geometry switch (0,1).\n");

    return_obj = make_unique<C3D4>(tag, std::move(node_tag), material_tag, is_true(nonlinear));
}

void new_c3d8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_c3d8() needs a valid tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_c3d8() needs eight valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_c3d8() needs a valid material tag.\n");
        return;
    }

    string reduced_scheme = "I";
    if(command.eof()) suanpan_debug("new_c3d8() assumes standard integration scheme (2*2).\n");
    else if(!get_input(command, reduced_scheme)) suanpan_error("new_c3d8() needs a valid reduced integration switch (0,1).\n");

    string nonlinear = "false";
    if(command.eof()) suanpan_debug("new_c3d8() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_c3d8() needs a valid nonlinear geometry switch (0,1).\n");

    return_obj = make_unique<C3D8>(tag, std::move(node_tag), material_tag, suanpan::to_upper(reduced_scheme[0]), is_true(nonlinear));
}

void new_c3d8r(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_c3d8r() needs a valid tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_c3d8r() needs eight valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_c3d8r() needs a valid material tag.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof()) suanpan_debug("new_c3d8() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_c3d8() needs a valid nonlinear geometry switch (0,1).\n");

    return_obj = make_unique<C3D8>(tag, std::move(node_tag), material_tag, 'R', is_true(nonlinear));
}

void new_c3d8i(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_c3d8i() needs a valid tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_c3d8i() needs eight valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_c3d8i() needs a valid material tag.\n");
        return;
    }

    return_obj = make_unique<C3D8I>(tag, std::move(node_tag), material_tag);
}

void new_cax3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cax3() needs a tag.\n");
        return;
    }

    uvec node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cax3() needs three valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cax3() needs a valid material tag.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof()) suanpan_debug("new_cax3() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_cax3() needs a valid nonlinear geometry switch.\n");

    return_obj = make_unique<CAX3>(tag, std::move(node_tag), material_tag, is_true(nonlinear));
}

void new_cax4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cax4() needs a tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cax4() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cax4() needs a valid material tag.\n");
        return;
    }

    return_obj = make_unique<CAX4>(tag, std::move(node_tag), material_tag, false);
}

void new_cax8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cax8() needs a tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cax8() needs eight valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cax8() needs a valid material tag.\n");
        return;
    }

    return_obj = make_unique<CAX8>(tag, std::move(node_tag), material_tag, false);
}

void new_contact2d(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_contact2d() needs a tag.\n");
        return;
    }

    unsigned master_tag, slave_tag;
    if(!get_input(command, master_tag)) {
        suanpan_error("new_contact2d() needs valid master group tag.\n");
        return;
    }
    if(!get_input(command, slave_tag)) {
        suanpan_error("new_contact2d() needs valid slave group tag.\n");
        return;
    }

    auto alpha = 1E6;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("new_contact2d() needs valid multiplier.\n");
        return;
    }

    return_obj = make_unique<Contact2D>(tag, master_tag, slave_tag, alpha);
}

void new_contact3d(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_contact3d() needs a tag.\n");
        return;
    }

    unsigned master_tag, slave_tag;
    if(!get_input(command, master_tag)) {
        suanpan_error("new_contact3d() needs valid master group tag.\n");
        return;
    }
    if(!get_input(command, slave_tag)) {
        suanpan_error("new_contact3d() needs valid slave group tag.\n");
        return;
    }

    auto alpha = 1E6;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("new_contact3d() needs valid multiplier.\n");
        return;
    }

    return_obj = make_unique<Contact3D>(tag, master_tag, slave_tag, alpha);
}

void new_cp3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cp3() needs a tag.\n");
        return;
    }

    uvec node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cp3() needs three valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cp3() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cp3() assumes thickness to be unit.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof()) suanpan_debug("new_cp3() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_cp3() needs a valid nonlinear geometry switch.\n");

    return_obj = make_unique<CP3>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
}

void new_cp4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cp4() needs a tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cp4() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cp4() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cp4() needs a valid thickness.\n");
        return;
    }

    string reduced_scheme = "N";
    if(command.eof()) suanpan_debug("new_cp4() assumes full integration.\n");
    else if(!get_input(command, reduced_scheme)) {
        suanpan_error("new_cp4() needs a valid reduced scheme switch.\n");
        return;
    }

    string nonlinear = "N";
    if(command.eof()) suanpan_debug("new_cp4() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) {
        suanpan_error("new_cp4() needs a valid nonlinear geometry switch.\n");
        return;
    }

    return_obj = make_unique<CP4>(tag, std::move(node_tag), material_tag, thickness, is_true(reduced_scheme), is_true(nonlinear));
}

void new_cp4i(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cp4i() needs a tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cp4i() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cp4i() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cp4i() needs a valid thickness.\n");
        return;
    }

    return_obj = make_unique<CP4I>(tag, std::move(node_tag), material_tag, thickness);
}

void new_cp4r(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cp4r() needs a tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cp4r() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cp4r() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cp4r() needs a valid thickness.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof()) suanpan_debug("new_cp4r() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_cp4r() needs a valid nonlinear geometry switch (0,1).\n");

    return_obj = make_unique<CP4>(tag, std::move(node_tag), material_tag, thickness, true, is_true(nonlinear));
}

void new_cp5(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cp5() needs a tag.\n");
        return;
    }

    uvec node_tag(5);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cp5() needs five valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cp5() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cp5() needs a valid thickness.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("new_cp5() needs a valid nonlinear geometry switch (0,1).\n");
        return;
    }

    return_obj = make_unique<CP5>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
}

void new_cp6(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cp6() needs a tag.\n");
        return;
    }

    uvec node_tag(6);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cp6() needs six valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cp6() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cp6() needs a valid thickness.\n");
        return;
    }

    string nonlinear = "false";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("new_cp6() needs a valid nonlinear geometry switch.\n");
        return;
    }

    return_obj = make_unique<CP6>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
}

void new_cp7(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cp7() needs a tag.\n");
        return;
    }

    uvec node_tag(7);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cp7() needs five valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cp7() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cp7() needs a valid thickness.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("new_cp7() needs a valid nonlinear geometry switch (0,1).\n");
        return;
    }

    return_obj = make_unique<CP7>(tag, std::move(node_tag), material_tag, thickness, is_true(nonlinear));
}

void new_cp8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cp8() needs a tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cp8() needs eight valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cp8() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cp8() needs a valid thickness.\n");
        return;
    }

    string reduced_scheme = "N";
    if(!get_optional_input(command, reduced_scheme)) {
        suanpan_error("new_cp8() needs a valid reduced integration switch.\n");
        return;
    }

    string nonlinear = "N";
    if(!get_optional_input(command, nonlinear)) {
        suanpan_error("new_cp8() needs a valid nonlinear geometry switch (0,1).\n");
        return;
    }

    return_obj = make_unique<CP8>(tag, std::move(node_tag), material_tag, thickness, is_true(reduced_scheme), is_true(nonlinear));
}

void new_cpe8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cpe8() needs a tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cpe8() needs eight valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cpe8() needs a valid material tag.\n");
        return;
    }

    string reduced_scheme = "N";
    if(!command.eof() && !get_input(command, reduced_scheme)) {
        suanpan_error("new_cp8() needs a valid reduced integration switch (0,1).\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("new_cp8() needs a valid nonlinear geometry switch (0,1).\n");
        return;
    }

    return_obj = make_unique<CP8>(tag, std::move(node_tag), material_tag, 1., is_true(reduced_scheme), is_true(nonlinear));
}

void new_cpe8r(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cpe8r() needs a tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cpe8r() needs eight valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cpe8r() needs a valid material tag.\n");
        return;
    }

    string nonlinear = "N";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("new_cpe8r() needs a valid nonlinear geometry switch.\n");
        return;
    }

    return_obj = make_unique<CP8>(tag, std::move(node_tag), material_tag, 1., true, is_true(nonlinear));
}

void new_cinp4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cinp4() needs a tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cinp4() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cinp4() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_cinp4() needs a valid thickness.\n");
        return;
    }

    return_obj = make_unique<CINP4>(tag, std::move(node_tag), material_tag, thickness);
}

void new_cin3d8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cin3d8() needs a tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_cin3d8() needs eight valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_cin3d8() needs a valid material tag.\n");
        return;
    }

    return_obj = make_unique<CIN3D8>(tag, std::move(node_tag), material_tag);
}

void new_csmt3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_csmt3() needs a tag.\n");
        return;
    }

    uvec node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_csmt3() needs three valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_csmt3() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_csmt3() needs a valid thickness.\n");
        return;
    }

    auto length = -1.;
    if(!get_optional_input(command, length)) {
        suanpan_error("new_csmt3() needs a valid length.\n");
        return;
    }

    return_obj = make_unique<CSMT3>(tag, std::move(node_tag), material_tag, thickness, length);
}

void new_csmt6(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_csmt6() needs a tag.\n");
        return;
    }

    uvec node_tag(6);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_csmt6() needs six valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_csmt6() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_csmt6() needs a valid thickness.\n");
        return;
    }

    auto length = -1.;
    if(!get_optional_input(command, length)) {
        suanpan_error("new_csmt6() needs a valid length.\n");
        return;
    }

    return_obj = make_unique<CSMT6>(tag, std::move(node_tag), material_tag, thickness, length);
}

void new_csmq(unique_ptr<Element>& return_obj, istringstream& command, const unsigned size) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_csmq() needs a tag.\n");
        return;
    }

    uvec node_tag(size);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_csmq() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_csmq() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_csmq() needs a valid thickness.\n");
        return;
    }

    auto length = -1.;
    if(!get_optional_input(command, length)) {
        suanpan_error("new_csmq() needs a valid length.\n");
        return;
    }

    if(4 == size) return_obj = make_unique<CSMQ4>(tag, std::move(node_tag), material_tag, thickness, length);
    else if(5 == size) return_obj = make_unique<CSMQ5>(tag, std::move(node_tag), material_tag, thickness, length);
    else if(6 == size) return_obj = make_unique<CSMQ6>(tag, std::move(node_tag), material_tag, thickness, length);
    else if(7 == size) return_obj = make_unique<CSMQ7>(tag, std::move(node_tag), material_tag, thickness, length);
    else if(8 == size) return_obj = make_unique<CSMQ8>(tag, std::move(node_tag), material_tag, thickness, length);
}

void new_damper01(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_damper01() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_damper01() needs two valid nodes.\n");
            return;
        }

    unsigned damper_tag;
    if(!get_input(command, damper_tag)) {
        suanpan_error("new_damper01() needs a valid tag.\n");
        return;
    }

    return_obj = make_unique<Damper01>(tag, std::move(node_tag), damper_tag);
}

void new_damper02(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_damper02() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_damper02() needs two valid nodes.\n");
            return;
        }

    unsigned damper_tag;
    if(!get_input(command, damper_tag)) {
        suanpan_error("new_damper02() needs a valid tag.\n");
        return;
    }

    unsigned spring_tag;
    if(!get_input(command, spring_tag)) {
        suanpan_error("new_damper02() needs a valid tag.\n");
        return;
    }

    string use_matrix = "true";
    if(!command.eof() && !get_input(command, use_matrix)) {
        suanpan_error("new_damper02() needs a valid switch.\n");
        return;
    }

    unsigned proceed = 0;
    if(!command.eof() && !get_input(command, proceed)) {
        suanpan_error("new_damper02() needs a valid proceed switch.\n");
        return;
    }

    auto beta = .5;
    if(!command.eof() && !get_input(command, beta)) {
        suanpan_error("new_damper02() needs a valid beta value.\n");
        return;
    }

    return_obj = make_unique<Damper02>(tag, std::move(node_tag), damper_tag, spring_tag, is_true(use_matrix), proceed, beta);
}

void new_dc3d4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_dc3d4() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_dc3d4() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_dc3d4() needs a valid material tag.\n");
        return;
    }

    double length;
    if(!get_input(command, length)) {
        suanpan_error("new_dc3d4() needs a valid characteristic length.\n");
        return;
    }

    double rate;
    if(!get_input(command, rate)) {
        suanpan_error("new_dc3d4() needs a valid energy release rate.\n");
        return;
    }

    return_obj = make_unique<DC3D4>(tag, std::move(node_tag), material_tag, length, rate);
}

void new_dc3d8(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_dc3d8() needs a valid tag.\n");
        return;
    }

    uvec node_tag(8);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_dc3d8() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_dc3d8() needs a valid material tag.\n");
        return;
    }

    double length;
    if(!get_input(command, length)) {
        suanpan_error("new_dc3d8() needs a valid characteristic length.\n");
        return;
    }

    double rate;
    if(!get_input(command, rate)) {
        suanpan_error("new_dc3d8() needs a valid energy release rate.\n");
        return;
    }

    return_obj = make_unique<DC3D8>(tag, std::move(node_tag), material_tag, length, rate);
}

void new_dcp3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_dcp3() needs a tag.\n");
        return;
    }

    uvec node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_dcp3() needs three valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_dcp3() needs a valid material tag.\n");
        return;
    }

    double length;
    if(!get_input(command, length)) {
        suanpan_error("new_dcp3() needs a valid characteristic length.\n");
        return;
    }

    double rate;
    if(!get_input(command, rate)) {
        suanpan_error("new_dcp3() needs a valid energy release rate.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_dcp3() needs a valid thickness.\n");
        return;
    }

    return_obj = make_unique<DCP3>(tag, std::move(node_tag), material_tag, length, rate, thickness);
}

void new_dcp4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_dcp4() needs a tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_dcp4() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_dcp4() needs a valid material tag.\n");
        return;
    }

    double length;
    if(!get_input(command, length)) {
        suanpan_error("new_dcp4() needs a valid characteristic length.\n");
        return;
    }

    double rate;
    if(!get_input(command, rate)) {
        suanpan_error("new_dcp4() needs a valid energy release rate.\n");
        return;
    }

    auto thickness = 1.;
    if(!get_optional_input(command, thickness)) {
        suanpan_error("new_dcp4() needs a valid thickness.\n");
        return;
    }

    return_obj = make_unique<DCP4>(tag, std::move(node_tag), material_tag, length, rate, thickness);
}

void new_dkt3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_dkt3() needs a valid tag.\n");
        return;
    }

    uvec node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_dkt3() needs three valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_dkt3() needs a valid material tag.\n");
        return;
    }

    double thickness;
    if(!get_input(command, thickness)) {
        suanpan_error("new_dkt3() needs a valid thickness.\n");
        return;
    }

    unsigned num_ip = 3;
    if(!get_optional_input(command, num_ip)) {
        suanpan_error("new_dkt3() needs a valid number of integration points.\n");
        return;
    }

    return_obj = make_unique<DKT3>(tag, std::move(node_tag), material_tag, thickness, num_ip);
}

void new_dkt4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_dkt4() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_dkt4() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_dkt4() needs a valid material tag.\n");
        return;
    }

    double thickness;
    if(!get_input(command, thickness)) {
        suanpan_error("new_dkt4() needs a valid thickness.\n");
        return;
    }

    unsigned num_ip = 3;
    if(!get_optional_input(command, num_ip)) {
        suanpan_error("new_dkt4() needs a valid number of integration points.\n");
        return;
    }

    return_obj = make_unique<DKT4>(tag, std::move(node_tag), material_tag, thickness, num_ip);
}

void new_dkts3(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_dkts3() needs a valid tag.\n");
        return;
    }

    uvec node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_dkts3() needs three valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_dkts3() needs a valid material tag.\n");
        return;
    }

    double thickness;
    if(!get_input(command, thickness)) {
        suanpan_error("new_dkts3() needs a valid thickness.\n");
        return;
    }

    unsigned num_ip = 3;
    if(!get_optional_input(command, num_ip)) {
        suanpan_error("new_dkts3() needs a valid number of integration points.\n");
        return;
    }

    return_obj = make_unique<DKTS3>(tag, std::move(node_tag), material_tag, thickness, num_ip);
}

void new_embedded(unique_ptr<Element>& return_obj, istringstream& command, const unsigned dof) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_embedded() needs a valid tag.\n");
        return;
    }

    unsigned element_tag;
    if(!get_input(command, element_tag)) {
        suanpan_error("new_embedded() needs a valid element tag.\n");
        return;
    }

    unsigned node_tag;
    if(!get_input(command, node_tag)) {
        suanpan_error("new_embedded() needs a valid node tag.\n");
        return;
    }

    auto alpha = 1E4;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("new_embedded() needs a valid alpha.\n");
        return;
    }

    if(2 == dof) return_obj = make_unique<Embedded2D>(tag, element_tag, node_tag, alpha);
    else return_obj = make_unique<Embedded3D>(tag, element_tag, node_tag, alpha);
}

void new_eb21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_eb21() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_eb21() needs two valid nodes.\n");
            return;
        }

    double area;
    if(!get_input(command, area)) {
        suanpan_error("new_eb21() needs a valid area.\n");
        return;
    }

    double moment_inertia;
    if(!get_input(command, moment_inertia)) {
        suanpan_error("new_eb21() needs a valid moment of inertia.\n");
        return;
    }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_eb21() needs a valid material tag.\n");
        return;
    }

    string nonlinear = "false";
    if(command.eof()) suanpan_debug("new_eb21() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_eb21() needs a valid nonlinear geometry switch (0,1).\n");

    return_obj = make_unique<EB21>(tag, std::move(node_tag), area, moment_inertia, material_tag, is_true(nonlinear));
}

void new_f21(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_f21() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_f21() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_f21() needs a valid section tag.\n");
        return;
    }

    unsigned int_pt = 6;
    if(!get_optional_input(command, int_pt)) {
        suanpan_error("new_f21() needs a valid number of integration points.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(command.eof()) suanpan_debug("new_f21() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_f21() needs a valid nonlinear geometry switch (0,1).\n");

    return_obj = make_unique<F21>(tag, std::move(node_tag), section_id, int_pt, !!nonlinear);
}

void new_f21h(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_f21h() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_f21h() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_f21h() needs a valid section tag.\n");
        return;
    }

    auto elastic_length = .2;
    if(!command.eof() && !get_input(command, elastic_length)) {
        suanpan_error("new_f21h() needs a valid number of integration points.\n");
        return;
    }

    unsigned nonlinear = 0;
    if(command.eof()) suanpan_debug("new_f21h() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_f21h() needs a valid nonlinear geometry switch (0,1).\n");

    return_obj = make_unique<F21H>(tag, std::move(node_tag), section_id, elastic_length, !!nonlinear);
}

void new_f31(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_f31() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_f31() needs two valid nodes.\n");
            return;
        }

    unsigned section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("new_f31() needs a valid section tag.\n");
        return;
    }

    unsigned orientation_id;
    if(!get_input(command, orientation_id)) {
        suanpan_error("new_f31() needs a valid orientation tag.\n");
        return;
    }

    unsigned int_pt = 5;
    if(!command.eof() && !get_input(command, int_pt)) {
        suanpan_error("new_f31() needs a valid number of integration points.\n");
        return;
    }

    string nonlinear = "false";
    if(!command.eof() && !get_input(command, nonlinear)) {
        suanpan_error("new_f31() needs a valid nonlinear geometry switch.\n");
        return;
    }

    return_obj = make_unique<F31>(tag, std::move(node_tag), section_id, orientation_id, int_pt, is_true(nonlinear));
}

void new_gcmq(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_gcmq() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_gcmq() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_gcmq() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_gcmq() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_gcmq() needs a valid thickness.\n");

    string int_scheme = "I";
    if(!command.eof() && !get_input(command, int_scheme)) suanpan_error("new_gcmq() needs a valid reduced scheme switch.\n");

    return_obj = make_unique<GCMQ>(tag, std::move(node_tag), material_tag, thickness, suanpan::to_upper(int_scheme[0]));
}

void new_gcmqi(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_gcmq() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_gcmq() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_gcmq() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_gcmq() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_gcmq() needs a valid thickness.\n");

    return_obj = make_unique<GCMQ>(tag, std::move(node_tag), material_tag, thickness, 'I');
}

void new_gcmql(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_gcmq() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_gcmq() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_gcmq() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_gcmq() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_gcmq() needs a valid thickness.\n");

    return_obj = make_unique<GCMQ>(tag, std::move(node_tag), material_tag, thickness, 'L');
}

void new_gcmqg(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_gcmq() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_gcmq() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_gcmq() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_gcmq() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_gcmq() needs a valid thickness.\n");

    return_obj = make_unique<GCMQ>(tag, std::move(node_tag), material_tag, thickness, 'G');
}

void new_sgcmqi(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_sgcmqi() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_sgcmqi() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_sgcmqi() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_sgcmqi() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_sgcmqi() needs a valid thickness.\n");

    return_obj = make_unique<SGCMQ>(tag, std::move(node_tag), material_tag, thickness, 'I');
}

void new_sgcmql(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_sgcmql() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_sgcmql() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_sgcmql() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_sgcmql() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_sgcmql() needs a valid thickness.\n");

    return_obj = make_unique<SGCMQ>(tag, std::move(node_tag), material_tag, thickness, 'L');
}

void new_sgcmqg(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_sgcmqg() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_sgcmqg() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_sgcmqg() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_sgcmqg() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_sgcmqg() needs a valid thickness.\n");

    return_obj = make_unique<SGCMQ>(tag, std::move(node_tag), material_tag, thickness, 'G');
}

void new_sgcms(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_sgcms() needs a tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_sgcms() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_sgcms() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && !get_input(command, thickness)) {
        suanpan_error("new_sgcms() needs a valid thickness.\n");
        return;
    }

    return_obj = make_unique<SGCMS>(tag, std::move(node_tag), material_tag, thickness);
}

void new_gq12(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_gq12() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_gq12() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_gq12() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_gq12() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_gq12() needs a valid thickness.\n");

    return_obj = make_unique<GQ12>(tag, std::move(node_tag), material_tag, thickness);
}

void new_joint(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_joint() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_joint() needs two valid nodes.\n");
            return;
        }

    vector<uword> pool;
    uword material_tag;
    while(!command.eof() && get_input(command, material_tag)) pool.emplace_back(material_tag);

    return_obj = make_unique<Joint>(tag, std::move(node_tag), std::move(pool));
}

void new_mass(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_mass() needs a valid tag.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("new_mass() needs one valid node.\n");
        return;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("new_mass() needs a valid magnitude.\n");
        return;
    }

    unsigned dof;
    vector<uword> dof_tag;
    while(!command.eof() && get_input(command, dof)) dof_tag.push_back(dof);

    return_obj = make_unique<Mass>(tag, node, magnitude, uvec(dof_tag));
}

void new_mindlin(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_mindlin() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_mindlin() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_mindlin() needs a valid material tag.\n");
        return;
    }

    double thickness;
    if(!get_input(command, thickness)) {
        suanpan_error("new_mindlin() needs a valid thickness.\n");
        return;
    }

    unsigned num_ip = 5;
    if(!command.eof() && !get_input(command, num_ip)) {
        suanpan_error("new_mindlin() needs a valid number of integration points.\n");
        return;
    }

    return_obj = make_unique<Mindlin>(tag, std::move(node_tag), material_tag, thickness, num_ip);
}

void new_mvlem(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_mvlem() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_mvlem() needs two valid nodes.\n");
            return;
        }

    unsigned shear_tag;
    if(!get_input(command, shear_tag)) {
        suanpan_error("new_mvlem() needs a valid shear material.\n");
        return;
    }

    double c_height;
    if(!get_input(command, c_height)) {
        suanpan_error("new_mvlem() needs a valid c.\n");
        return;
    }

    vector<double> B, H, R;
    vector<uword> CT, ST;
    while(!command.eof()) {
        double t_value;
        uword t_tag;
        if(!get_input(command, t_value)) {
            suanpan_error("new_mvlem() needs a valid fibre width.\n");
            return;
        }
        B.emplace_back(t_value);
        if(!get_input(command, t_value)) {
            suanpan_error("new_mvlem() needs a valid fibre thickness.\n");
            return;
        }
        H.emplace_back(t_value);
        if(!get_input(command, t_value)) {
            suanpan_error("new_mvlem() needs a valid fibre reinforcement ratio.\n");
            return;
        }
        R.emplace_back(t_value);
        if(!get_input(command, t_tag)) {
            suanpan_error("new_mvlem() needs a valid material tag.\n");
            return;
        }
        CT.emplace_back(t_tag);
        if(!get_input(command, t_tag)) {
            suanpan_error("new_mvlem() needs a valid material tag.\n");
            return;
        }
        ST.emplace_back(t_tag);
    }

    return_obj = make_unique<MVLEM>(tag, std::move(node_tag), B, H, R, uvec(CT), uvec(ST), shear_tag, c_height);
}

void new_pcpedc(unique_ptr<Element>& return_obj, istringstream& command, const unsigned node) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_pcpedc() needs a tag.\n");
        return;
    }

    uvec node_tag(node);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_pcpedc() needs %u valid nodes.\n", node);
            return;
        }

    unsigned solid_tag;
    if(!get_input(command, solid_tag)) {
        suanpan_error("new_pcpedc() needs a valid material tag for solid phase.\n");
        return;
    }
    unsigned fluid_tag;
    if(!get_input(command, fluid_tag)) {
        suanpan_error("new_pcpedc() needs a valid material tag for fluid phase.\n");
        return;
    }

    double alpha, n, k;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("new_pcpedc() needs a valid alpha.\n");
        return;
    }
    if(!get_optional_input(command, n)) {
        suanpan_error("new_pcpedc() needs a valid porosity.\n");
        return;
    }
    if(!get_optional_input(command, k)) {
        suanpan_error("new_pcpedc() needs a valid permeability.\n");
        return;
    }

    if(4 == node) return_obj = make_unique<PCPE4DC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n, k);
    else if(8 == node) return_obj = make_unique<PCPE8DC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n, k);
}

void new_pcpeuc(unique_ptr<Element>& return_obj, istringstream& command, const unsigned node) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_pcpeuc() needs a tag.\n");
        return;
    }

    uvec node_tag(node);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_pcpeuc() needs %u valid nodes.\n", node);
            return;
        }

    unsigned solid_tag;
    if(!get_input(command, solid_tag)) {
        suanpan_error("new_pcpeuc() needs a valid material tag for solid phase.\n");
        return;
    }
    unsigned fluid_tag;
    if(!get_input(command, fluid_tag)) {
        suanpan_error("new_pcpeuc() needs a valid material tag for fluid phase.\n");
        return;
    }

    double alpha, n;
    if(!get_optional_input(command, alpha)) {
        suanpan_error("new_pcpeuc() needs a valid alpha.\n");
        return;
    }
    if(!get_optional_input(command, n)) {
        suanpan_error("new_pcpeuc() needs a valid porosity.\n");
        return;
    }

    if(4 == node) return_obj = make_unique<PCPE4UC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n);
    else if(8 == node) return_obj = make_unique<PCPE8UC>(tag, std::move(node_tag), solid_tag, fluid_tag, alpha, n);
}

void new_ps(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_ps() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_ps() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_ps() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_ps() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_ps() needs a valid thickness.\n");

    return_obj = make_unique<PS>(tag, std::move(node_tag), material_tag, thickness);
}

void new_qe2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_qe2() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_qe2() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_qe2() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(command.eof()) suanpan_debug("new_qe2() assumes thickness to be unit.\n");
    else if(!get_input(command, thickness)) suanpan_error("new_qe2() needs a valid thickness.\n");

    return_obj = make_unique<QE2>(tag, std::move(node_tag), material_tag, thickness);
}

void new_s4(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_s4() needs a tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_s4() needs four valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_s4() needs a valid material tag.\n");
        return;
    }

    auto thickness = 1.;
    if(!command.eof() && !get_input(command, thickness)) {
        suanpan_error("new_s4() needs a valid thickness.\n");
        return;
    }

    return_obj = make_unique<S4>(tag, std::move(node_tag), material_tag, thickness);
}

void new_singlesection2d(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_singlesection2d() needs a valid tag.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("new_singlesection2d() needs one valid node.\n");
        return;
    }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("new_singlesection2d() needs a valid section tag.\n");
        return;
    }

    return_obj = make_unique<SingleSection2D>(tag, node, section_tag);
}

void new_singlesection3d(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_singlesection3d() needs a valid tag.\n");
        return;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("new_singlesection3d() needs one valid node.\n");
        return;
    }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("new_singlesection3d() needs a valid section tag.\n");
        return;
    }

    return_obj = make_unique<SingleSection3D>(tag, node, section_tag);
}

void new_spring01(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_spring01() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_spring01() needs two valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_spring01() needs a valid material tag.\n");
        return;
    }

    return_obj = make_unique<Spring01>(tag, std::move(node_tag), material_tag);
}

void new_spring02(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_spring02() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_spring02() needs two valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_spring02() needs a valid material tag.\n");
        return;
    }

    return_obj = make_unique<Spring02>(tag, std::move(node_tag), material_tag);
}

void new_t2d2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_t2d2() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_t2d2() needs two valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_t2d2() needs a valid material tag.\n");
        return;
    }

    double area;
    if(!get_input(command, area)) {
        suanpan_error("new_t2d2() needs a valid area.\n");
        return;
    }

    string nonlinear = "N", update_area = "N", log_strain = "N";

    if(command.eof()) suanpan_debug("new_t2d2() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_t2d2() needs a valid nonlinear geometry switch.\n");

    if(command.eof()) suanpan_debug("new_t2d2() assumes constant area.\n");
    else if(!get_input(command, update_area)) suanpan_error("new_t2d2() needs a valid switch (0,1) to indicate if to update area.\n");

    if(command.eof()) suanpan_debug("new_t2d2() assumes engineering strain.\n");
    else if(!get_input(command, log_strain)) suanpan_error("new_t2d2() needs a valid switch (0,1) to indicate if to use engineering strain.\n");

    return_obj = make_unique<T2D2>(tag, std::move(node_tag), material_tag, area, is_true(nonlinear), is_true(update_area), is_true(log_strain));
}

void new_t2d2s(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_t2d2s() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_t2d2s() needs two valid nodes.\n");
            return;
        }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("new_t2d2s() needs a valid material/section tag.\n");
        return;
    }

    string nonlinear = "N", log_strain = "N";

    if(command.eof()) suanpan_debug("new_t2d2() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_t2d2s() needs a valid nonlinear geometry switch.\n");

    if(command.eof()) suanpan_debug("new_t2d2() assumes engineering strain.\n");
    else if(!get_input(command, log_strain)) suanpan_error("new_t2d2s() needs a valid switch to indicate if to use engineering strain.\n");

    return_obj = make_unique<T2D2S>(tag, std::move(node_tag), section_tag, is_true(nonlinear), is_true(log_strain));
}

void new_t3d2(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_t3d2() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_t3d2() needs two valid nodes.\n");
            return;
        }

    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("new_t3d2() needs a valid material tag.\n");
        return;
    }

    double area;
    if(!get_input(command, area)) {
        suanpan_error("new_t3d2() needs a valid area.\n");
        return;
    }

    string nonlinear = "N", update_area = "N", log_strain = "N";

    if(command.eof()) suanpan_debug("new_t2d2() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_t3d2() needs a valid nonlinear geometry switch (0,1).\n");

    if(command.eof()) suanpan_debug("new_truss2d() assumes constant area.\n");
    else if(!get_input(command, update_area)) suanpan_error("new_t3d2() needs a valid switch (0,1) to indicate if update area.\n");

    if(command.eof()) suanpan_debug("new_t3d2() assumes engineering strain.\n");
    else if(!get_input(command, log_strain)) suanpan_error("new_t3d2() needs a valid switch (0,1) to indicate if to use engineering strain.\n");

    return_obj = make_unique<T3D2>(tag, std::move(node_tag), material_tag, area, is_true(nonlinear), is_true(update_area), is_true(log_strain));
}

void new_t3d2s(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_t3d2s() needs a valid tag.\n");
        return;
    }

    uvec node_tag(2);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_t3d2s() needs two valid nodes.\n");
            return;
        }

    unsigned section_tag;
    if(!get_input(command, section_tag)) {
        suanpan_error("new_t3d2s() needs a valid material/section tag.\n");
        return;
    }

    string nonlinear = "N", log_strain = "N";

    if(command.eof()) suanpan_debug("new_t2d2() assumes linear geometry.\n");
    else if(!get_input(command, nonlinear)) suanpan_error("new_t3d2s() needs a valid nonlinear geometry switch.\n");

    if(command.eof()) suanpan_debug("new_t2d2() assumes engineering strain.\n");
    else if(!get_input(command, log_strain)) suanpan_error("new_t3d2s() needs a valid switch to indicate if to use engineering strain.\n");

    return_obj = make_unique<T3D2S>(tag, std::move(node_tag), section_tag, is_true(nonlinear), is_true(log_strain));
}

void new_tie(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_tie() needs a valid tag.\n");
        return;
    }

    double magnitude, penalty;
    if(!get_input(command, magnitude) || !get_input(command, penalty)) {
        suanpan_error("new_tie() needs a valid magnitude.\n");
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

void new_patchquad(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_patchquad() needs a tag.\n");
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
                    suanpan_error("new_patchquad() needs a valid material tag.\n");
                    return;
                }
                ignore_whitespace(command);
            }
            else if(is_equal(parameter, "-thickness")) {
                if(!get_input(command, thickness)) {
                    suanpan_error("new_patchquad() needs a valid thickness.\n");
                    return;
                }
                ignore_whitespace(command);
            }
        }
    }

    if(node_tag.empty()) {
        suanpan_error("new_patchquad() needs a valid node tags.\n");
        return;
    }
    if(knot_x.empty()) {
        suanpan_error("new_patchquad() needs a valid knot vector.\n");
        return;
    }
    if(knot_y.empty()) {
        suanpan_error("new_patchquad() needs a valid knot vector.\n");
        return;
    }
    if(-1 == material_tag) {
        suanpan_error("new_patchquad() needs a valid material tag.\n");
        return;
    }

    return_obj = make_unique<PatchQuad>(tag, std::move(knot_x), std::move(knot_y), std::move(node_tag), static_cast<unsigned>(material_tag), thickness);
}

void new_patchcube(unique_ptr<Element>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_patchcube() needs a tag.\n");
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
                    suanpan_error("new_patchcube() needs a valid material tag.\n");
                    return;
                }
                ignore_whitespace(command);
            }
        }
    }

    if(node_tag.empty()) {
        suanpan_error("new_patchcube() needs a valid node tags.\n");
        return;
    }
    if(knot_x.empty()) {
        suanpan_error("new_patchcube() needs a valid knot vector.\n");
        return;
    }
    if(knot_y.empty()) {
        suanpan_error("new_patchcube() needs a valid knot vector.\n");
        return;
    }
    if(knot_z.empty()) {
        suanpan_error("new_patchcube() needs a valid knot vector.\n");
        return;
    }
    if(-1 == material_tag) {
        suanpan_error("new_patchcube() needs a valid material tag.\n");
        return;
    }

    return_obj = make_unique<PatchCube>(tag, std::move(knot_x), std::move(knot_y), std::move(knot_z), std::move(node_tag), static_cast<unsigned>(material_tag));
}
