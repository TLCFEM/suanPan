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

// ReSharper disable StringLiteralTypo
// ReSharper disable IdentifierTypo
#include "LoadParser.h"

#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Load/Load>
#include <Toolbox/resampling.h>

namespace {
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

        unsigned dof_id;
        if(!get_input(command, dof_id)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        return_obj = std::make_unique<NodalAcceleration>(load_id, magnitude, get_remaining<uword>(command), dof_id, amplitude_id);
    }

    void new_bodyforce(unique_ptr<Load>& return_obj, std::istringstream& command, const bool flag) {
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

        unsigned dof_id;
        if(!get_input(command, dof_id)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(flag) return_obj = std::make_unique<GroupBodyForce>(load_id, magnitude, get_remaining<uword>(command), uvec{dof_id}, amplitude_id);
        else return_obj = std::make_unique<BodyForce>(load_id, magnitude, get_remaining<uword>(command), uvec{dof_id}, amplitude_id);
    }

    void new_cload(unique_ptr<Load>& return_obj, std::istringstream& command, const bool flag) {
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

        unsigned dof_id;
        if(!get_input(command, dof_id)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(flag) return_obj = std::make_unique<GroupNodalForce>(load_id, magnitude, get_remaining<uword>(command), dof_id, amplitude_id);
        else return_obj = std::make_unique<NodalForce>(load_id, magnitude, get_remaining<uword>(command), dof_id, amplitude_id);
    }

    void new_refload(unique_ptr<Load>& return_obj, std::istringstream& command) {
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

        unsigned dof_id;
        if(!get_input(command, dof_id)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        return_obj = std::make_unique<ReferenceForce>(load_id, magnitude, get_remaining<uword>(command), dof_id);
    }

    void new_lineudl(unique_ptr<Load>& return_obj, std::istringstream& command, const unsigned dimension) {
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

        unsigned dof_id;
        if(!get_input(command, dof_id)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(2 == dimension) return_obj = std::make_unique<LineUDL2D>(load_id, magnitude, get_remaining<uword>(command), dof_id, amplitude_id);
        else return_obj = std::make_unique<LineUDL3D>(load_id, magnitude, get_remaining<uword>(command), dof_id, amplitude_id);
    }

    void new_displacement(unique_ptr<Load>& return_obj, std::istringstream& command, const bool flag) {
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

        unsigned dof_id;
        if(!get_input(command, dof_id)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(flag) return_obj = std::make_unique<GroupNodalDisplacement>(load_id, magnitude, get_remaining<uword>(command), uvec{dof_id}, amplitude_id);
        else return_obj = std::make_unique<NodalDisplacement>(load_id, magnitude, get_remaining<uword>(command), uvec{dof_id}, amplitude_id);
    }

    void new_supportmotion(unique_ptr<Load>& return_obj, std::istringstream& command, const unsigned flag) {
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

        unsigned dof_id;
        if(!get_input(command, dof_id)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        if(0 == flag) return_obj = std::make_unique<SupportDisplacement>(load_id, magnitude, get_remaining<uword>(command), dof_id, amplitude_id);
        else if(1 == flag) return_obj = std::make_unique<SupportVelocity>(load_id, magnitude, get_remaining<uword>(command), dof_id, amplitude_id);
        else return_obj = std::make_unique<SupportAcceleration>(load_id, magnitude, get_remaining<uword>(command), dof_id, amplitude_id);
    }
} // namespace

int create_new_amplitude(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string amplitude_type;
    if(!get_input(command, amplitude_type)) {
        suanpan_error("A valid amplitude type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid amplitude type is required.\n");
        return SUANPAN_SUCCESS;
    }

    if(const auto step_tag = domain->get_current_step_tag(); is_equal(amplitude_type, "Constant")) domain->insert(std::make_shared<Constant>(tag, step_tag));
    else if(is_equal(amplitude_type, "Ramp")) domain->insert(std::make_shared<Ramp>(tag, step_tag));
    else if(is_equal(amplitude_type, "Tabular")) {
        std::string file_name;
        if(!get_input(command, file_name)) {
            suanpan_error("A valid file is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(command.eof()) domain->insert(std::make_shared<Tabular>(tag, std::move(file_name), step_tag));
        else {
            uword up_rate;
            if(!get_input(command, up_rate)) {
                suanpan_error("A valid upsampling ratio is required.\n");
                return SUANPAN_SUCCESS;
            }

            std::string window_type = "Hamming";
            if(!get_optional_input(command, window_type)) {
                suanpan_error("A valid window type is required.\n");
                return SUANPAN_SUCCESS;
            }

            auto window_size = 8llu;
            if(!get_optional_input(command, window_size)) {
                suanpan_error("A valid window size is required.\n");
                return SUANPAN_SUCCESS;
            }

            const auto result = upsampling(window_type, file_name, up_rate, window_size);

            if(result.empty()) {
                suanpan_error("Fail to perform upsampling.\n");
                return SUANPAN_SUCCESS;
            }

            domain->insert(std::make_shared<Tabular>(tag, result.col(0), result.col(1), step_tag));
        }
    }
    else if(is_equal(amplitude_type, "TabularSpline")) {
        std::string file_name;
        if(!get_input(command, file_name)) {
            suanpan_error("A valid file is required.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(std::make_shared<TabularSpline>(tag, std::move(file_name), step_tag));
    }
    else if(is_equal(amplitude_type, "Decay")) {
        double A, TD;
        if(!get_input(command, A, TD)) {
            suanpan_error("A valid value is required.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(std::make_shared<Decay>(tag, A, TD, step_tag));
    }
    else if(is_equal(amplitude_type, "Linear")) {
        double A;
        if(!get_input(command, A)) {
            suanpan_error("A valid value is required.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(std::make_shared<Linear>(tag, A, step_tag));
    }
    else if(is_equal(amplitude_type, "Combine")) {
        domain->insert(std::make_shared<Combine>(tag, get_remaining<uword>(command), step_tag));
    }
    else if(is_equal(amplitude_type, "Custom")) {
        unsigned expression;
        if(!get_input(command, expression)) {
            suanpan_error("A valid expression tag is required.\n");
            return SUANPAN_SUCCESS;
        }
        domain->insert(std::make_shared<CustomAmplitude>(tag, expression, step_tag));
    }
    else if(is_equal(amplitude_type, "Modulated") || is_equal(amplitude_type, "Sine") || is_equal(amplitude_type, "Cosine")) {
        double W;
        if(!get_input(command, W)) {
            suanpan_error("A valid value is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(amplitude_type, "Modulated")) domain->insert(std::make_shared<Modulated>(tag, W, get_remaining<double>(command), step_tag));
        else if(is_equal(amplitude_type, "Sine")) domain->insert(std::make_shared<Sine>(tag, W, get_remaining<double>(command), step_tag));
        else if(is_equal(amplitude_type, "Cosine")) domain->insert(std::make_shared<Cosine>(tag, W, get_remaining<double>(command), step_tag));
    }
    else if(is_equal(amplitude_type, "NZStrongMotion")) {
        std::string name;
        if(!get_input(command, name)) {
            suanpan_error("A valid name is required.\n");
            return SUANPAN_SUCCESS;
        }

        domain->insert(std::make_shared<NZStrongMotion>(tag, name.c_str(), step_tag));
    }

    return SUANPAN_SUCCESS;
}

int create_new_load(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string load_id;
    if(!get_input(command, load_id)) {
        suanpan_error("A valid load type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Load> new_load = nullptr;

    if(is_equal(load_id, "Acceleration")) new_acceleration(new_load, command);
    else if(is_equal(load_id, "BodyForce")) new_bodyforce(new_load, command, false);
    else if(is_equal(load_id, "Cload")) new_cload(new_load, command, false);
    else if(is_equal(load_id, "Disp") || is_equal(load_id, "Displacement") || is_equal(load_id, "DispLoad")) new_displacement(new_load, command, false);
    else if(is_equal(load_id, "GroupBodyForce")) new_bodyforce(new_load, command, true);
    else if(is_equal(load_id, "GroupCload")) new_cload(new_load, command, true);
    else if(is_equal(load_id, "GroupDisp") || is_equal(load_id, "GroupDisplacement") || is_equal(load_id, "GroupDispLoad")) new_displacement(new_load, command, true);
    else if(is_equal(load_id, "LineUDL2D")) new_lineudl(new_load, command, 2);
    else if(is_equal(load_id, "LineUDL3D")) new_lineudl(new_load, command, 3);
    else if(is_equal(load_id, "ReferenceLoad") || is_equal(load_id, "RefLoad") || is_equal(load_id, "RefForce")) new_refload(new_load, command);
    else if(is_equal(load_id, "SupportAcceleration")) new_supportmotion(new_load, command, 2);
    else if(is_equal(load_id, "SupportDisplacement")) new_supportmotion(new_load, command, 0);
    else if(is_equal(load_id, "SupportVelocity")) new_supportmotion(new_load, command, 1);
    else external_module::object(new_load, domain, load_id, command);

    if(new_load != nullptr) new_load->set_start_step(domain->get_current_step_tag());

    if(new_load == nullptr || !domain->insert(std::move(new_load)))
        suanpan_error("Fail to create new load via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
