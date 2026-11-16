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

#include "RecorderParser.h"

#include <Domain/ConditionalModifier.h>
#include <Domain/DomainBase.h>
#include <Recorder/Recorder>
#include <Toolbox/utility.h>

namespace {
    struct RecorderArguments {
        OutputType variable_token{OutputType::NL};
        unsigned interval{1u};
    };

    bool get_recorder_arguments(std::istringstream& command, RecorderArguments& arguments, const bool with_variable = true) {
        if(with_variable) {
            std::string variable_type;
            if(!get_input(command, variable_type)) {
                suanpan_error("A valid recorder type is required.\n");
                return false;
            }
            arguments.variable_token = to_token(std::move(variable_type));
        }

        while(true)
            if(const auto peek_value = command.peek(); is_equal(peek_value, '\t') || is_equal(peek_value, ' ')) command.ignore();
            else break;

        if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
            std::string tmp_string;
            get_input(command, tmp_string);
            if(!get_input(command, arguments.interval)) return false;
        }

        return true;
    }

    void new_eigen_recorder(unique_ptr<Recorder>& return_obj, const unsigned tag, const bool use_hdf5) {
        return_obj = std::make_unique<EigenRecorder>(tag, use_hdf5);
    }

    void new_modal_recorder(unique_ptr<Recorder>& return_obj, std::istringstream& command, const unsigned tag, const bool use_hdf5) {
        std::vector<Node::DOF> dof_pool;
        for(const auto& token : get_remaining<std::string>(command)) {
            const auto dof = parse_dof(token);
            if(dof.empty()) {
                suanpan_error("A valid dof identifier is required.\n");
                return;
            }
            dof_pool.insert(dof_pool.end(), dof.cbegin(), dof.cend());
        }

        if(dof_pool.empty()) {
            suanpan_error("At least one dof identifier is required.\n");
            return;
        }

        return_obj = std::make_unique<ModalRecorder>(tag, std::move(dof_pool), use_hdf5);
    }

    void new_modal_force_recorder(unique_ptr<Recorder>& return_obj, const unsigned tag, const bool use_hdf5) {
        return_obj = std::make_unique<ModalForceRecorder>(tag, use_hdf5);
    }

    void new_frame_recorder(unique_ptr<Recorder>& return_obj, std::istringstream& command, const unsigned tag) {
        RecorderArguments arguments;
        if(!get_recorder_arguments(command, arguments)) return;

        return_obj = std::make_unique<FrameRecorder>(tag, arguments.variable_token, arguments.interval);
    }

    void new_visualisation_recorder(unique_ptr<Recorder>& return_obj, std::istringstream& command, const unsigned tag) {
        RecorderArguments arguments;
        if(!get_recorder_arguments(command, arguments)) return;

        std::string para;
        auto width = 6u;
        auto scale = 1.;
        while(!command.eof() && get_input(command, para))
            if(is_equal(para, "Width")) {
                if(!get_input(command, width)) suanpan_error("A valid width is required, using six.\n");
            }
            else if(is_equal(para, "Scale")) {
                if(!get_input(command, scale)) suanpan_error("A valid scale is required, using unity.\n");
            }

        return_obj = std::make_unique<VisualisationRecorder>(tag, arguments.variable_token, arguments.interval, width, scale);
#ifndef SUANPAN_VTK
        suanpan_warning("Visualisation recorders have no effect as the required VTK support is not available.\n");
#endif
    }

    template<typename T> void new_tag_recorder(unique_ptr<Recorder>& return_obj, std::istringstream& command, const unsigned tag, const bool use_hdf5) {
        RecorderArguments arguments;
        if(!get_recorder_arguments(command, arguments)) return;

        return_obj = std::make_unique<T>(tag, uvec(get_remaining<uword>(command)), arguments.variable_token, arguments.interval, use_hdf5);
    }

    void new_amplitude_recorder(unique_ptr<Recorder>& return_obj, std::istringstream& command, const unsigned tag, const bool use_hdf5) {
        RecorderArguments arguments;
        if(!get_recorder_arguments(command, arguments, false)) return;

        return_obj = std::make_unique<AmplitudeRecorder>(tag, uvec(get_remaining<uword>(command)), OutputType::AMP, arguments.interval, use_hdf5);
    }

    void new_global_recorder(unique_ptr<Recorder>& return_obj, std::istringstream& command, const unsigned tag, const bool use_hdf5) {
        RecorderArguments arguments;
        if(!get_recorder_arguments(command, arguments)) return;

        if(OutputType::K == arguments.variable_token) return_obj = std::make_unique<GlobalStiffnessRecorder>(tag, arguments.interval, use_hdf5);
        else if(OutputType::M == arguments.variable_token) return_obj = std::make_unique<GlobalMassRecorder>(tag, arguments.interval, use_hdf5);
        else return_obj = std::make_unique<GlobalRecorder>(tag, arguments.variable_token, arguments.interval, use_hdf5);
    }

    int dispatch_new_recorder(const shared_ptr<DomainBase>& domain, std::istringstream& command, const unsigned tag, const bool use_hdf5) {
        std::string object_type;
        if(!get_input(command, object_type)) {
            suanpan_error("A valid object type is required.\n");
            return SUANPAN_SUCCESS;
        }

        unique_ptr<Recorder> new_recorder = nullptr;

        if(is_equal(object_type, "Eigen")) new_eigen_recorder(new_recorder, tag, use_hdf5);
        else if(is_equal(object_type, "Modal")) new_modal_recorder(new_recorder, command, tag, use_hdf5);
        else if(is_equal(object_type, "ModalForce")) new_modal_force_recorder(new_recorder, tag, use_hdf5);
        else if(is_equal(object_type, "Frame")) new_frame_recorder(new_recorder, command, tag);
        else if(is_equal(object_type, "Visualisation")) new_visualisation_recorder(new_recorder, command, tag);
        else if(is_equal(object_type, "Node")) new_tag_recorder<NodeRecorder>(new_recorder, command, tag, use_hdf5);
        else if(is_equal(object_type, "GroupNode")) new_tag_recorder<GroupNodeRecorder>(new_recorder, command, tag, use_hdf5);
        else if(is_equal(object_type, "Sum")) new_tag_recorder<SumRecorder>(new_recorder, command, tag, use_hdf5);
        else if(is_equal(object_type, "GroupSum")) new_tag_recorder<GroupSumRecorder>(new_recorder, command, tag, use_hdf5);
        else if(is_equal(object_type, "Element")) new_tag_recorder<ElementRecorder>(new_recorder, command, tag, use_hdf5);
        else if(is_equal(object_type, "GroupElement")) new_tag_recorder<GroupElementRecorder>(new_recorder, command, tag, use_hdf5);
        else if(is_equal(object_type, "Amplitude")) new_amplitude_recorder(new_recorder, command, tag, use_hdf5);
        else if(is_equal(object_type, "Global")) new_global_recorder(new_recorder, command, tag, use_hdf5);

        if(nullptr == new_recorder || !domain->insert(std::move(new_recorder))) suanpan_error("Fail to create new recorder via \"{}\".\n", command.str());

        return SUANPAN_SUCCESS;
    }
} // namespace

int create_new_recorder(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::string file_type;
    if(!get_input(command, file_type)) {
        suanpan_error("A valid object type is required.\n");
        return SUANPAN_SUCCESS;
    }

    return dispatch_new_recorder(domain, command, tag, if_startswith(file_type, "hdf"));
}

int create_new_recorder(const shared_ptr<DomainBase>& domain, std::istringstream& command, const bool use_hdf5) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    return dispatch_new_recorder(domain, command, tag, use_hdf5);
}
