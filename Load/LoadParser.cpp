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
#include "LoadParser.h"

#include <Domain/ExternalModule.h>
#include <Load/Load>
#include <Toolbox/resampling.h>

namespace {
    //
    // loads
    //

    void new_acceleration(unique_ptr<Load>& return_obj, std::istringstream& command) {
        unsigned load_id;
        if(!get_input(command, load_id)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned amplitude_id;
        if(!get_input(command, amplitude_id)) {
            suanpan_error("A valid amplitude tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid load magnitude is required.\n");
            return;
        }

        std::string dof_token;
        if(!get_input(command, dof_token)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        auto dof_pool = parse_dof(dof_token);
        if(dof_pool.empty()) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        return_obj = std::make_unique<NodalAcceleration>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
    }

    template<bool flag> void new_bodyforce(unique_ptr<Load>& return_obj, std::istringstream& command) {
        unsigned load_id;
        if(!get_input(command, load_id)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned amplitude_id;
        if(!get_input(command, amplitude_id)) {
            suanpan_error("A valid amplitude tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid load magnitude is required.\n");
            return;
        }

        std::string dof_token;
        if(!get_input(command, dof_token)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        auto dof_pool = parse_dof(dof_token);
        if(dof_pool.empty()) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(flag) return_obj = std::make_unique<GroupBodyForce>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
        else return_obj = std::make_unique<BodyForce>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
    }

    template<bool flag> void new_cload(unique_ptr<Load>& return_obj, std::istringstream& command) {
        unsigned load_id;
        if(!get_input(command, load_id)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned amplitude_id;
        if(!get_input(command, amplitude_id)) {
            suanpan_error("A valid amplitude tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid load magnitude is required.\n");
            return;
        }

        std::string dof_token;
        if(!get_input(command, dof_token)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        auto dof_pool = parse_dof(dof_token);
        if(dof_pool.empty()) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(flag) return_obj = std::make_unique<GroupNodalForce>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
        else return_obj = std::make_unique<NodalForce>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
    }

    template<bool flag> void new_refload(unique_ptr<Load>& return_obj, std::istringstream& command) {
        unsigned load_id;
        if(!get_input(command, load_id)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        if(unsigned amplitude_id; !get_input(command, amplitude_id)) {
            suanpan_error("A valid amplitude tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid load magnitude is required.\n");
            return;
        }

        std::string dof_token;
        if(!get_input(command, dof_token)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        auto dof_pool = parse_dof(dof_token);
        if(dof_pool.empty()) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(flag) return_obj = std::make_unique<GroupReferenceForce>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool));
        else return_obj = std::make_unique<ReferenceForce>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool));
    }

    template<unsigned dimension> void new_lineudl(unique_ptr<Load>& return_obj, std::istringstream& command) {
        unsigned load_id;
        if(!get_input(command, load_id)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned amplitude_id;
        if(!get_input(command, amplitude_id)) {
            suanpan_error("A valid amplitude tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid load magnitude is required.\n");
            return;
        }

        std::string dof_token;
        if(!get_input(command, dof_token)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        auto dof_pool = parse_dof(dof_token);
        if(1 != dof_pool.size()) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(2 == dimension) return_obj = std::make_unique<LineUDL2D>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
        else return_obj = std::make_unique<LineUDL3D>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
    }

    template<bool flag> void new_displacement(unique_ptr<Load>& return_obj, std::istringstream& command) {
        unsigned load_id;
        if(!get_input(command, load_id)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned amplitude_id;
        if(!get_input(command, amplitude_id)) {
            suanpan_error("A valid amplitude tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid load magnitude is required.\n");
            return;
        }

        std::string dof_token;
        if(!get_input(command, dof_token)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        auto dof_pool = parse_dof(dof_token);
        if(dof_pool.empty()) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(flag) return_obj = std::make_unique<GroupNodalDisplacement>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
        else return_obj = std::make_unique<NodalDisplacement>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
    }

    template<unsigned flag> void new_supportmotion(unique_ptr<Load>& return_obj, std::istringstream& command) {
        unsigned load_id;
        if(!get_input(command, load_id)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned amplitude_id;
        if(!get_input(command, amplitude_id)) {
            suanpan_error("A valid amplitude tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid load magnitude is required.\n");
            return;
        }

        std::string dof_token;
        if(!get_input(command, dof_token)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        auto dof_pool = parse_dof(dof_token);
        if(dof_pool.empty()) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(0 == flag) return_obj = std::make_unique<SupportDisplacement>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
        else if(1 == flag) return_obj = std::make_unique<SupportVelocity>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
        else return_obj = std::make_unique<SupportAcceleration>(load_id, magnitude, get_remaining<uword>(command), std::move(dof_pool), amplitude_id);
    }

    //
    // amplitudes
    //

    template<typename T> void new_tag_amplitude(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag);
    }

    void new_tabular_amplitude(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        std::string file_name;
        if(!get_input(command, file_name)) {
            suanpan_error("A valid file is required.\n");
            return;
        }

        if(command.eof()) {
            return_obj = std::make_unique<Tabular>(tag, std::move(file_name));
            return;
        }

        uword up_rate;
        if(!get_input(command, up_rate)) {
            suanpan_error("A valid upsampling ratio is required.\n");
            return;
        }

        std::string window_type = "Hamming";
        if(!get_optional_input(command, window_type)) {
            suanpan_error("A valid window type is required.\n");
            return;
        }

        uword window_size = 8u;
        if(!get_optional_input(command, window_size)) {
            suanpan_error("A valid window size is required.\n");
            return;
        }

        const auto result = upsampling(window_type, file_name, up_rate, window_size);

        if(result.empty()) {
            suanpan_error("Fail to perform upsampling.\n");
            return;
        }

        return_obj = std::make_unique<Tabular>(tag, result.col(0), result.col(1));
    }

    void new_tabular_spline_amplitude(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        std::string file_name;
        if(!get_input(command, file_name)) {
            suanpan_error("A valid file is required.\n");
            return;
        }

        return_obj = std::make_unique<TabularSpline>(tag, std::move(file_name));
    }

    void new_decay_amplitude(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double A, TD;
        if(!get_input(command, A, TD)) {
            suanpan_error("A valid value is required.\n");
            return;
        }

        return_obj = std::make_unique<Decay>(tag, A, TD);
    }

    void new_linear_amplitude(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double A;
        if(!get_input(command, A)) {
            suanpan_error("A valid value is required.\n");
            return;
        }

        return_obj = std::make_unique<Linear>(tag, A);
    }

    void new_combine_amplitude(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Combine>(tag, get_remaining<uword>(command));
    }

    void new_custom_amplitude(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned expression;
        if(!get_input(command, expression)) {
            suanpan_error("A valid expression tag is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomAmplitude>(tag, expression);
    }

    template<typename T> void new_trigonometric_amplitude(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double W;
        if(!get_input(command, W)) {
            suanpan_error("A valid value is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag, W, get_remaining<double>(command));
    }

    void new_nz_strong_motion(unique_ptr<Amplitude>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        std::string name;
        if(!get_input(command, name)) {
            suanpan_error("A valid name is required.\n");
            return;
        }

        return_obj = std::make_unique<NZStrongMotion>(tag, name.c_str());
    }
} // namespace

int create_new_load(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("A valid load type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Load> new_load = nullptr;

    if(is_equal(load_id, "Acceleration")) new_acceleration(new_load, command);
    else if(is_equal(load_id, "BodyForce")) new_bodyforce<false>(new_load, command);
    else if(is_equal(load_id, "Cload")) new_cload<false>(new_load, command);
    else if(is_equal_any(load_id, "Disp", "Displacement", "DispLoad")) new_displacement<false>(new_load, command);
    else if(is_equal(load_id, "GroupBodyForce")) new_bodyforce<true>(new_load, command);
    else if(is_equal(load_id, "GroupCload")) new_cload<true>(new_load, command);
    else if(is_equal_any(load_id, "GroupDisp", "GroupDisplacement", "GroupDispLoad")) new_displacement<true>(new_load, command);
    else if(is_equal(load_id, "LineUDL2D")) new_lineudl<2>(new_load, command);
    else if(is_equal(load_id, "LineUDL3D")) new_lineudl<3>(new_load, command);
    else if(is_equal_any(load_id, "ReferenceLoad", "RefLoad", "RefForce")) new_refload<false>(new_load, command);
    else if(is_equal_any(load_id, "GroupReferenceLoad", "GroupRefLoad", "GroupRefForce")) new_refload<true>(new_load, command);
    else if(is_equal(load_id, "SupportAcceleration")) new_supportmotion<2>(new_load, command);
    else if(is_equal(load_id, "SupportDisplacement")) new_supportmotion<0>(new_load, command);
    else if(is_equal(load_id, "SupportVelocity")) new_supportmotion<1>(new_load, command);
    else external_module::object(new_load, domain, load_id, command);

    if(new_load != nullptr) new_load->set_start_step(domain->get_current_step_tag());

    if(new_load == nullptr || !domain->insert(std::move(new_load)))
        suanpan_error("Fail to create new load via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_amplitude(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string amplitude_type;
    if(!get_input(command, amplitude_type)) {
        suanpan_error("A valid amplitude type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Amplitude> new_amplitude = nullptr;

    if(is_equal(amplitude_type, "Constant")) new_tag_amplitude<Constant>(new_amplitude, command);
    else if(is_equal(amplitude_type, "Ramp")) new_tag_amplitude<Ramp>(new_amplitude, command);
    else if(is_equal(amplitude_type, "Tabular")) new_tabular_amplitude(new_amplitude, command);
    else if(is_equal(amplitude_type, "TabularSpline")) new_tabular_spline_amplitude(new_amplitude, command);
    else if(is_equal(amplitude_type, "Decay")) new_decay_amplitude(new_amplitude, command);
    else if(is_equal(amplitude_type, "Linear")) new_linear_amplitude(new_amplitude, command);
    else if(is_equal(amplitude_type, "Combine")) new_combine_amplitude(new_amplitude, command);
    else if(is_equal(amplitude_type, "Custom")) new_custom_amplitude(new_amplitude, command);
    else if(is_equal(amplitude_type, "Modulated")) new_trigonometric_amplitude<Modulated>(new_amplitude, command);
    else if(is_equal(amplitude_type, "Sine")) new_trigonometric_amplitude<Sine>(new_amplitude, command);
    else if(is_equal(amplitude_type, "Cosine")) new_trigonometric_amplitude<Cosine>(new_amplitude, command);
    else if(is_equal(amplitude_type, "NZStrongMotion")) new_nz_strong_motion(new_amplitude, command);

    if(nullptr == new_amplitude || !domain->insert(std::move(new_amplitude))) suanpan_error("Fail to create new amplitude via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
