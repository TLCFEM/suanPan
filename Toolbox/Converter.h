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

#ifndef CONVERTER_H
#define CONVERTER_H

#include <fstream>
#include <unordered_map>

struct Converter {
    static std::istringstream clean(std::string);
    static std::string extract_name(const std::string&, const std::string& = "");

    int boundary_tag = 0;
    int set_tag = 0;
    int amplitude_tag = 0;
    int load_tag = 0;
    int step_tag = 0;
    int integrator_tag = 0;
    int converger_tag = 0;

    std::string output_path;

    std::unordered_map<std::string, int> node_set_pool;
    std::unordered_map<std::string, int> element_set_pool;
    std::unordered_map<std::string, int> amplitude_pool;
    std::unordered_map<std::string, int> step_pool;

    explicit Converter(std::string&&);

    int process(std::ifstream&, std::ofstream&);

    int process_amplitude(std::ifstream&, std::ofstream&, const std::string&, const std::string&);
    int process_boundary(std::ifstream&, std::ofstream&, const std::string&, const std::string&);
    int process_cload(std::ifstream&, std::ofstream&, const std::string&);
    int process_element(std::ifstream&, std::ofstream&, const std::string&) const;
    int process_node(std::ifstream&, std::ofstream&) const;
    int process_node_set(std::ifstream&, std::ofstream&, const std::string&, bool);
    int process_element_set(std::ifstream&, std::ofstream&, const std::string&, bool);
    int process_step(std::ifstream&, std::ofstream&, const std::string&);

    [[nodiscard]] std::string extract_element_type(const std::string&) const;

    bool getline(std::istream&, std::string&) const;
};

#endif
