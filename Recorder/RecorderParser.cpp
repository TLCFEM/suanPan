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

#include "RecorderParser.h"
#include <Domain/DomainBase.h>
#include <Recorder/Recorder>
#include <Toolbox/utility.h>

int create_new_recorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        SP_E("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    string file_type;
    if(!get_input(command, file_type)) {
        SP_E("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        SP_E("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Eigen")) {
        if(!domain->insert(make_shared<EigenRecorder>(tag, is_equal(file_type[0], 'h')))) SP_E("Fail to create new eigen recorder.\n");
        return SUANPAN_SUCCESS;
    }

    string variable_type;
    if(!is_equal(object_type, "Amplitude") && !get_input(command, variable_type)) {
        SP_E("A valid recorder type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned interval = 1;

    while(true)
        if(const auto peek_value = command.peek(); is_equal(peek_value, '\t') || is_equal(peek_value, ' ')) command.ignore();
        else break;

    if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
        string tmp_string;
        get_input(command, tmp_string);
        if(!get_input(command, interval)) return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Frame")) {
        if(!domain->insert(make_shared<FrameRecorder>(tag, to_list(variable_type.c_str()), interval))) SP_E("Fail to create new frame recorder.\n");
        return SUANPAN_SUCCESS;
    }
    if(is_equal(object_type, "Visualisation")) {
        unsigned width = 6;
        if(!command.eof() && !get_input(command, width)) width = 6;
        if(!domain->insert(make_shared<VisualisationRecorder>(tag, to_list(variable_type.c_str()), interval, width))) SP_E("Fail to create new visualisation recorder.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned s_object_tag;
    std::vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(const auto use_hdf5 = is_equal(file_type[0], 'h'); is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5)))
        SP_E("Fail to create new node recorder.\n");
    else if(is_equal(object_type, "GroupNode") && !domain->insert(make_shared<GroupNodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5)))
        SP_E("Fail to create new group node recorder.\n");
    else if(is_equal(object_type, "Sum") && !domain->insert(make_shared<SumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5)))
        SP_E("Fail to create new summation recorder.\n");
    else if(is_equal(object_type, "GroupSum") && !domain->insert(make_shared<GroupSumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5)))
        SP_E("Fail to create new group summation recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5)))
        SP_E("Fail to create new element recorder.\n");
    else if(is_equal(object_type, "GroupElement") && !domain->insert(make_shared<GroupElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, use_hdf5)))
        SP_E("Fail to create new group element recorder.\n");
    else if(is_equal(object_type, "Amplitude") && !domain->insert(make_shared<AmplitudeRecorder>(tag, uvec(object_tag), OutputType::AMP, interval, true, use_hdf5)))
        SP_E("Fail to create new amplitude recorder.\n");
    else if(is_equal(object_type, "Global")) {
        bool flag;
        if(OutputType::K == to_list(variable_type.c_str()))
            flag = domain->insert(make_shared<GlobalStiffnessRecorder>(tag, interval, true, use_hdf5));
        else if(OutputType::M == to_list(variable_type.c_str()))
            flag = domain->insert(make_shared<GlobalMassRecorder>(tag, interval, true, use_hdf5));
        else
            flag = domain->insert(make_shared<GlobalRecorder>(tag, to_list(variable_type.c_str()), interval, true, use_hdf5));
        if(!flag) SP_E("Fail to create new global recorder.\n");
    }

    return SUANPAN_SUCCESS;
}

int create_new_plainrecorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        SP_E("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        SP_E("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Eigen")) {
        if(!domain->insert(make_shared<EigenRecorder>(tag, false))) SP_E("Fail to create new eigen recorder.\n");
        return SUANPAN_SUCCESS;
    }

    string variable_type;
    if(!is_equal(object_type, "Amplitude") && !get_input(command, variable_type)) {
        SP_E("A valid recorder type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned interval = 1;

    while(true)
        if(const auto peek_value = command.peek(); is_equal(peek_value, '\t') || is_equal(peek_value, ' ')) command.ignore();
        else break;

    if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
        string tmp_string;
        get_input(command, tmp_string);
        if(!get_input(command, interval)) return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Visualisation")) {
        unsigned width = 6;
        if(!command.eof() && !get_input(command, width)) width = 6;
        if(!domain->insert(make_shared<VisualisationRecorder>(tag, to_list(variable_type.c_str()), interval, width))) SP_E("Fail to create new visualisation recorder.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned s_object_tag;
    std::vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false)))
        SP_E("Fail to create new node recorder.\n");
    else if(is_equal(object_type, "GroupNode") && !domain->insert(make_shared<GroupNodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false)))
        SP_E("Fail to create new group node recorder.\n");
    else if(is_equal(object_type, "Sum") && !domain->insert(make_shared<SumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false)))
        SP_E("Fail to create new summation recorder.\n");
    else if(is_equal(object_type, "GroupSum") && !domain->insert(make_shared<GroupSumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false)))
        SP_E("Fail to create new group summation recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false)))
        SP_E("Fail to create new element recorder.\n");
    else if(is_equal(object_type, "GroupElement") && !domain->insert(make_shared<GroupElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, false)))
        SP_E("Fail to create new group element recorder.\n");
    else if(is_equal(object_type, "Amplitude") && !domain->insert(make_shared<AmplitudeRecorder>(tag, uvec(object_tag), OutputType::AMP, interval, true, false)))
        SP_E("Fail to create new amplitude recorder.\n");
    else if(is_equal(object_type, "Global")) {
        bool flag;
        if(OutputType::K == to_list(variable_type.c_str()))
            flag = domain->insert(make_shared<GlobalStiffnessRecorder>(tag, interval, true, false));
        else if(OutputType::M == to_list(variable_type.c_str()))
            flag = domain->insert(make_shared<GlobalMassRecorder>(tag, interval, true, false));
        else
            flag = domain->insert(make_shared<GlobalRecorder>(tag, to_list(variable_type.c_str()), interval, true, false));
        if(!flag) SP_E("Fail to create new global recorder.\n");
    }

    return SUANPAN_SUCCESS;
}

int create_new_hdf5recorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        SP_E("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    string object_type;
    if(!get_input(command, object_type)) {
        SP_E("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Eigen")) {
        if(!domain->insert(make_shared<EigenRecorder>(tag, true))) SP_E("Fail to create new eigen recorder.\n");
        return SUANPAN_SUCCESS;
    }

    string variable_type;
    if(!is_equal(object_type, "Amplitude") && !get_input(command, variable_type)) {
        SP_E("A valid recorder type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned interval = 1;

    while(true)
        if(const auto peek_value = command.peek(); is_equal(peek_value, '\t') || is_equal(peek_value, ' ')) command.ignore();
        else break;

    if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
        string tmp_string;
        get_input(command, tmp_string);
        if(!get_input(command, interval)) return SUANPAN_SUCCESS;
    }

    if(is_equal(object_type, "Frame")) {
        if(!domain->insert(make_shared<FrameRecorder>(tag, to_list(variable_type.c_str()), interval))) SP_E("Fail to create new frame recorder.\n");
        return SUANPAN_SUCCESS;
    }
    if(is_equal(object_type, "Visualisation")) {
        string para;
        unsigned width = 6;
        auto scale = 1.;
        while(!command.eof() && get_input(command, para))
            if(is_equal(para, "Width")) {
                if(!get_input(command, width)) {
                    width = 6;
                    SP_E("A valid width is required.\n");
                }
            }
            else if(is_equal(para, "Scale")) {
                if(!get_input(command, scale)) {
                    scale = 1.;
                    SP_E("A valid scale is required.\n");
                }
            }
        if(!domain->insert(make_shared<VisualisationRecorder>(tag, to_list(variable_type.c_str()), interval, width, scale))) SP_E("Fail to create new visualisation recorder.\n");
        return SUANPAN_SUCCESS;
    }

    uword s_object_tag = 0;
    std::vector<uword> object_tag;
    while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

    if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true)))
        SP_E("Fail to create new node recorder.\n");
    else if(is_equal(object_type, "GroupNode") && !domain->insert(make_shared<GroupNodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true)))
        SP_E("Fail to create new group node recorder.\n");
    else if(is_equal(object_type, "Sum") && !domain->insert(make_shared<SumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true)))
        SP_E("Fail to create new summation recorder.\n");
    else if(is_equal(object_type, "GroupSum") && !domain->insert(make_shared<GroupSumRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true)))
        SP_E("Fail to create new group summation recorder.\n");
    else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true)))
        SP_E("Fail to create new element recorder.\n");
    else if(is_equal(object_type, "GroupElement") && !domain->insert(make_shared<GroupElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), interval, true, true)))
        SP_E("Fail to create new group element recorder.\n");
    else if(is_equal(object_type, "Amplitude") && !domain->insert(make_shared<AmplitudeRecorder>(tag, uvec(object_tag), OutputType::AMP, interval, true, true)))
        SP_E("Fail to create new amplitude recorder.\n");
    else if(is_equal(object_type, "Global")) {
        bool flag;
        if(OutputType::K == to_list(variable_type.c_str()))
            flag = domain->insert(make_shared<GlobalStiffnessRecorder>(tag, interval, true, true));
        else if(OutputType::M == to_list(variable_type.c_str()))
            flag = domain->insert(make_shared<GlobalMassRecorder>(tag, interval, true, true));
        else
            flag = domain->insert(make_shared<GlobalRecorder>(tag, to_list(variable_type.c_str()), interval, true, true));
        if(!flag) SP_E("Fail to create new global recorder.\n");
    }

    return SUANPAN_SUCCESS;
}
