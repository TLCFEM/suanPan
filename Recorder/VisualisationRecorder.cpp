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

void VisualisationRecorder::record_impl([[maybe_unused]] const shared_ptr<DomainBase>& D) {
#ifdef SUANPAN_VTK
    std::ostringstream file_name;

    file_name << 'R' << get_tag() << '-' << to_name(original_type) << '-' << std::setw(width) << std::setfill('0') << ++total_counter;

    const fs::path file_path = SUANPAN_OUTPUT;

    config.file_name = (file_path / file_name.str()).generic_string();

    vtk_cell_plot(D, config);
#endif
}

VisualisationRecorder::VisualisationRecorder(const unsigned T, const OutputType L, const unsigned I, const int W, [[maybe_unused]] const double S)
    : Recorder(T, {}, L, I, false)
    , width(W) {
#ifdef SUANPAN_VTK
    config.set(L);
    config.scale = S;
#endif
}

void VisualisationRecorder::clear_status() {
    total_counter = 0u;
    Recorder::clear_status();
}

void VisualisationRecorder::save() {}

void VisualisationRecorder::print() {
    suanpan_info("A visualisation recorder.\n");
}
