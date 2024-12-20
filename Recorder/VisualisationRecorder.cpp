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

#include "VisualisationRecorder.h"
#include <Domain/DomainBase.h>

extern fs::path SUANPAN_OUTPUT;

VisualisationRecorder::VisualisationRecorder(const unsigned T, const OutputType L, const unsigned I, const unsigned W, [[maybe_unused]] const double S)
    : Recorder(T, {}, L, I, false, false)
    , width(W) {
#ifdef SUANPAN_VTK
    config.save_file = true;
    config.type = get_variable_type();
    config.scale = S;

    const auto P = to_token(to_category(config.type));

    function_handler = OutputType::U == P || OutputType::V == P || OutputType::A == P || OutputType::RF == P || OutputType::DF == P || OutputType::IF == P ? &vtk_plot_node_quantity : &vtk_plot_element_quantity;
#endif
}

void VisualisationRecorder::record([[maybe_unused]] const shared_ptr<DomainBase>& D) {
#ifdef SUANPAN_VTK
    if(!if_perform_record()) return;

    ostringstream file_name;

    file_name << 'R' << get_tag() << '-' << to_name(get_variable_type()) << '-' << std::setw(static_cast<int>(width)) << std::setfill('0') << ++total_counter << ".vtk";

    fs::path file_path = SUANPAN_OUTPUT;

    file_path.append(file_name.str());

    config.file_name = file_path.generic_string();

    (*function_handler)(D, config);
#endif
}

void VisualisationRecorder::save() {}

void VisualisationRecorder::print() {
    suanpan_info("A visualisation recorder.\n");
}
