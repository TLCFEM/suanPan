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

#include "Converter.h"
#include <Toolbox/utility.h>
#include <iomanip>
#include <sstream>

Converter::Converter(std::string&& OP)
    : output_path(std::forward<string>(OP)) {}

int Converter::process(std::ifstream& input_file, std::ofstream& output_file) {
    output_file.setf(std::ios::scientific);
    output_file << std::setprecision(3);

    std::string line;

    while(getline(input_file, line)) {
        if(if_contain(line, "**")) continue;
        if(if_contain(line, "*node")) {
            //
            if(SUANPAN_SUCCESS != process_node(input_file, output_file)) break;
        }
        else if(if_contain(line, "*element")) {
            //
            if(SUANPAN_SUCCESS != process_element(input_file, output_file, extract_element_type(line))) break;
        }
        else if(if_contain(line, "*boundary")) {
            //
            if(SUANPAN_SUCCESS != process_boundary(input_file, output_file, extract_name(line, "amplitude"), extract_name(line, "type"))) break;
        }
        else if(if_contain(line, "*step")) {
            //
            if(SUANPAN_SUCCESS != process_step(input_file, output_file, extract_name(line, "name"))) break;
        }
        else if(if_contain(line, "*amplitude")) {
            //
            if(SUANPAN_SUCCESS != process_amplitude(input_file, output_file, extract_name(line, "name"), extract_name(line, "definition"))) break;
        }
        else if(if_contain(line, "*cload")) {
            //
            if(SUANPAN_SUCCESS != process_cload(input_file, output_file, extract_name(line, "amplitude"))) break;
        }
        else if(if_contain(line, "*nset") && if_contain(line, "instance")) {
            // only deal with sets created in instance
            if(SUANPAN_SUCCESS != process_node_set(input_file, output_file, extract_name(line, "nset"), if_contain(line, "generate"))) break;
        }
        else if(if_contain(line, "*elset") && if_contain(line, "instance")) {
            // only deal with sets created in instance
            if(SUANPAN_SUCCESS != process_element_set(input_file, output_file, extract_name(line, "elset"), if_contain(line, "generate"))) break;
        }
    }

    return SUANPAN_SUCCESS;
}

int Converter::process_amplitude(std::ifstream& input_file, std::ofstream& output_file, const std::string& amplitude_name, const std::string& type) {
    if(std::string line; is_equal(type, "")) {
        if(std::ofstream amp_out(output_path + amplitude_name); amp_out.is_open()) {
            amp_out.setf(std::ios::scientific);
            amp_out << std::setprecision(3);
            amplitude_pool.emplace(amplitude_name, ++amplitude_tag);
            output_file << "amplitude Tabular " << amplitude_tag << ' ' << amplitude_name << '\n';
            while(input_file.peek() != '*' && getline(input_file, line)) {
                double value;
                double time;
                auto tmp_str = clean(line);
                while(get_input(tmp_str, time) && get_input(tmp_str, value)) amp_out << time << '\t' << value << '\n';
            }
        }
    }
    else if(is_equal(type, "modulated")) {
        amplitude_pool.emplace(amplitude_name, ++amplitude_tag);
        output_file << "amplitude Modulated " << amplitude_tag << ' ';
        std::vector<double> pool;
        pool.reserve(5);
        if(getline(input_file, line)) {
            auto tmp_str = clean(line);
            for(auto I = 0; I < 5; ++I) {
                double value;
                get_input(tmp_str, value);
                pool.emplace_back(value);
            }
            output_file << pool[1] << ' ' << .5 * pool[3] / std::acos(-1.) << ' ' << .5 * pool[4] / std::acos(-1.) << '\n';
            suanpan_warning("The starting time is set to the start of current step.\n");
            if(pool[0] != 0.)
                suanpan_warning("The constant part is not supported currently.\n");
        }
    }
    else if(is_equal(type, "decay")) {
        amplitude_pool.emplace(amplitude_name, ++amplitude_tag);
        output_file << "amplitude Decay " << amplitude_tag << ' ';
        std::vector<double> pool;
        pool.reserve(4);
        if(getline(input_file, line)) {
            auto tmp_str = clean(line);
            for(auto I = 0; I < 4; ++I) {
                double value;
                get_input(tmp_str, value);
                pool.emplace_back(value);
            }
            output_file << pool[1] << ' ' << pool[3] << '\n';
            suanpan_warning("The starting time is set to the start of current step.\n");
            if(pool[0] != 0.)
                suanpan_warning("The constant part is not supported currently.\n");
        }
    }
    else if(is_equal(type, "periodic")) {
        suanpan_warning("The periodic amplitude is not supported currently.\n");
        return SUANPAN_SUCCESS;
    }

    output_file << '\n';

    return SUANPAN_SUCCESS;
}

int Converter::process_boundary(std::ifstream& input_file, std::ofstream& output_file, const std::string& amplitude_name, const std::string& type_name) {
    std::string command_name;

    while(input_file.peek() != '*')
        if(std::string line; getline(input_file, line)) {
            auto tmp_str = clean(line);
            std::vector<std::string> pool;
            pool.reserve(4);
            std::string item;
            while(get_input(tmp_str, item)) pool.emplace_back(item);
            auto t_set_tag = 0;
            if(pool.size() > 1) {
                if(node_set_pool.contains(pool[0])) {
                    // find corresponding node set
                    t_set_tag = node_set_pool[pool[0]];
                }
                else {
                    // directly apply on node
                    get_input(tmp_str = std::istringstream(pool[0]), t_set_tag);
                }
            }
            if(4 == pool.size()) {
                // node set, first dof, last dof, nonzero value ==> load
                int dof_a, dof_b;
                double value;
                get_input(tmp_str = std::istringstream(pool[1]), dof_a);
                get_input(tmp_str = std::istringstream(pool[2]), dof_b);
                get_input(tmp_str = std::istringstream(pool[3]), value);

                int amp_tag;
                if(amplitude_name.empty() || !amplitude_pool.contains(amplitude_name)) amp_tag = 0;
                else amp_tag = amplitude_pool[amplitude_name];

                if(type_name.empty() || if_contain(type_name, "displacement")) command_name = "groupdisplacement ";
                else if(if_contain(type_name, "velocity")) command_name = "groupvelocity ";
                else if(if_contain(type_name, "acceleration")) command_name = "groupacceleration ";

                while(dof_a <= dof_b) output_file << command_name << ++load_tag << ' ' << amp_tag << ' ' << value << ' ' << dof_a++ << ' ' << t_set_tag << '\n';
            }
            else if(3 == pool.size()) {
                if(if_contain(pool[2], ".")) {
                    // node set, dof, nonzero double >> load
                    int dof;
                    double value;
                    get_input(tmp_str = std::istringstream(pool[1]), dof);
                    get_input(tmp_str = std::istringstream(pool[2]), value);

                    int amp_tag;
                    if(amplitude_name.empty() || !amplitude_pool.contains(amplitude_name)) amp_tag = 0;
                    else amp_tag = amplitude_pool[amplitude_name];

                    if(type_name.empty() || if_contain(type_name, "displacement")) command_name = "groupdisplacement ";
                    else if(if_contain(type_name, "velocity")) command_name = "groupvelocity ";
                    else if(if_contain(type_name, "acceleration")) command_name = "groupacceleration ";

                    output_file << command_name << ++load_tag << ' ' << amp_tag << ' ' << value << ' ' << dof << ' ' << t_set_tag << '\n';
                }
                else {
                    // node set, first dof, last dof ==> boundary condition
                    int dof_a, dof_b;
                    get_input(tmp_str = std::istringstream(pool[1]), dof_a);
                    get_input(tmp_str = std::istringstream(pool[2]), dof_b);
                    while(dof_a <= dof_b) output_file << "groupmultiplierbc " << ++boundary_tag << ' ' << dof_a++ << ' ' << t_set_tag << '\n';
                }
            }
            else if(2 == pool.size()) {
                // node set, dof ==> boundary condition
                output_file << "groupmultiplierbc " << ++boundary_tag << ' ' << pool[1] << ' ' << t_set_tag << '\n';
            }
        }

    output_file << '\n';

    return SUANPAN_SUCCESS;
}

int Converter::process_cload(std::ifstream& input_file, std::ofstream& output_file, const std::string& amplitude_name) {
    while(input_file.peek() != '*')
        if(std::string line; getline(input_file, line)) {
            output_file << "groupcload " << ++load_tag << ' ' << (amplitude_name.empty() ? 0 : amplitude_pool[amplitude_name]) << ' ';
            auto t_str = clean(line);
            if(std::string r_set_tag, dof, value; get_input(t_str, r_set_tag) && get_input(t_str, dof) && get_input(t_str, value)) output_file << value << ' ' << dof << ' ' << node_set_pool[r_set_tag];
            output_file << '\n';
        }

    return SUANPAN_SUCCESS;
}

int Converter::process_element(std::ifstream& input_file, std::ofstream& output_file, const std::string& type) const {
    while(input_file.peek() != '*')
        if(std::string line; getline(input_file, line)) {
            unsigned tag;
            auto tmp_str = clean(line);
            output_file << "element " << type;
            while(get_input(tmp_str, tag)) output_file << ' ' << tag;
            output_file << '\n';
        }

    output_file << '\n';

    return SUANPAN_SUCCESS;
}

std::string Converter::extract_element_type(const std::string& line) const {
    auto element_type = extract_name(line, "type");

    if(is_equal(element_type, "CPS4R")) return "CP4R";
    if(is_equal(element_type, "CPS4")) return "CP4";
    if(is_equal(element_type, "CPS8")) return "CP8";
    if(is_equal(element_type, "CPS3")) return "CP3";
    if(is_equal(element_type, "CPS6")) return "CP6";
    if(is_equal(element_type, "CPE4R")) return "CP4R";
    if(is_equal(element_type, "CPE4")) return "CP4";
    if(is_equal(element_type, "CPE8")) return "CP8";
    if(is_equal(element_type, "CPE3")) return "CP3";
    if(is_equal(element_type, "CPE6")) return "CP6";
    if(is_equal(element_type, "C3D8")) return "C3D8";
    if(is_equal(element_type, "C3D8R")) return "C3D8R";

    return element_type;
}

int Converter::process_node(std::ifstream& input_file, std::ofstream& output_file) const {
    while(input_file.peek() != '*')
        if(std::string line; getline(input_file, line)) {
            double component;
            unsigned tag;
            auto tmp_str = clean(line);
            if(!get_input(tmp_str, tag)) return SUANPAN_SUCCESS;
            output_file << "node " << tag;
            while(get_input(tmp_str, component)) output_file << ' ' << component;
            output_file << '\n';
        }

    output_file << '\n';

    return SUANPAN_SUCCESS;
}

int Converter::process_node_set(std::ifstream& input_file, std::ofstream& output_file, const std::string& set_name, const bool generate) {
    while(input_file.peek() != '*')
        if(std::string line; getline(input_file, line)) {
            node_set_pool.emplace(set_name, ++set_tag);
            if(generate) output_file << "generate ";
            output_file << "nodegroup " << set_tag;
            if(generate) {
                auto tmp_str = clean(line);
                if(int start; get_input(tmp_str, start)) {
                    output_file << ' ' << start;
                    if(int end, interval; get_input(tmp_str, end) && get_input(tmp_str, interval)) output_file << ' ' << interval << ' ' << end;
                }
            }
            else
                do {
                    int tag;
                    auto tmp_str = clean(line);
                    while(get_input(tmp_str, tag)) output_file << ' ' << tag;
                }
                while(input_file.peek() != '*' && getline(input_file, line));
            output_file << '\n';
        }

    output_file << '\n';

    return SUANPAN_SUCCESS;
}

int Converter::process_element_set(std::ifstream& input_file, std::ofstream& output_file, const std::string& set_name, const bool generate) {
    while(input_file.peek() != '*')
        if(std::string line; getline(input_file, line)) {
            element_set_pool.emplace(set_name, ++set_tag);
            if(generate) output_file << "generate ";
            output_file << "elementgroup " << set_tag;
            if(generate) {
                auto tmp_str = clean(line);
                if(int start; get_input(tmp_str, start)) {
                    output_file << ' ' << start;
                    if(int end, interval; get_input(tmp_str, end) && get_input(tmp_str, interval)) output_file << ' ' << interval << ' ' << end;
                }
            }
            else
                do {
                    int tag;
                    auto tmp_str = clean(line);
                    while(get_input(tmp_str, tag)) output_file << ' ' << tag;
                }
                while(input_file.peek() != '*' && getline(input_file, line));
            output_file << '\n';
        }

    output_file << '\n';

    return SUANPAN_SUCCESS;
}

int Converter::process_step(std::ifstream& input_file, std::ofstream& output_file, const std::string& step_name) {
    if(std::string line; getline(input_file, line)) {
        if(auto tmp_str = clean(line); is_equal(line.substr(0, 7), "*Static")) {
            step_pool.emplace(step_name, ++step_tag);
            output_file << "step static " << step_tag;
            if(getline(input_file, line)) {
                tmp_str = clean(line);
                std::vector<double> pool;
                pool.reserve(5);
                double value;
                while(get_input(tmp_str, value)) pool.emplace_back(value);
                output_file << ' ' << pool[1] << "\nset ini_step_size " << pool[0] << '\n';
                if(2 == pool.size()) output_file << "set fixed_step_size true\n";
                else if(4 == pool.size()) output_file << "set min_step_size " << pool[2] << "\nset max_step_size " << pool[3] << '\n';
            }
        }
        else if(is_equal(line.substr(0, 8), "*Dynamic")) {
            step_pool.emplace(step_name, ++step_tag);
            output_file << "step dynamic " << step_tag;
            if(getline(input_file, line)) {
                tmp_str = clean(line);
                std::vector<double> pool;
                pool.reserve(5);
                double value;
                while(get_input(tmp_str, value)) pool.emplace_back(value);
                output_file << ' ' << pool[1] << "\nset ini_step_size " << pool[0] << '\n';
                if(2 == pool.size()) output_file << "set fixed_step_size true\n";
                else if(3 == pool.size()) output_file << "set min_step_size " << pool[2] << '\n';
            }

            output_file << "\nintegrator " << ++integrator_tag << " Newmark\n";
        }

        output_file << "\nconverger AbsIncreDisp " << ++converger_tag << " 1E-8 7 1\n";
    }

    output_file << '\n';

    return SUANPAN_SUCCESS;
}

bool Converter::getline(std::istream& in, std::string& out) const {
    do {
        std::getline(in, out);
        if(in.fail()) return false;
    }
    while(out.empty());

    for(auto& I : out) I = static_cast<char>(std::tolower(static_cast<int>(I)));

    return true;
}

std::istringstream Converter::clean(std::string line) {
    for(auto& c : line) if(',' == c) c = ' ';

    return std::istringstream(line);
}

std::string Converter::extract_name(const std::string& line, const std::string& key) {
    auto s_pos = line.find(key + '=');

    // in case of no name
    if(std::string::npos == s_pos) return "";

    auto e_pos = (s_pos += 1 + key.length());
    while(e_pos != line.length() && line[e_pos] != ',') ++e_pos;

    return line.substr(s_pos, e_pos - s_pos);
}
