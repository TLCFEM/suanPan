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

#include "FrameRecorder.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Element/Element.h>

extern fs::path SUANPAN_OUTPUT;

void FrameRecorder::record_impl([[maybe_unused]] const shared_ptr<DomainBase>& D) {
#ifdef SUANPAN_HDF5
    std::ostringstream group_name;
    group_name << '/' << D->get_factory()->get_current_time();

    auto group = file.create_group(group_name.str());

    for(auto& I : D->get_element_pool()) {
        if(const auto data = normalise_size(I->record(variable_type)); !data.empty()) {
            mat data_to_write(data[0].n_elem, static_cast<uword>(data.size()));

            uword idx = 0;
            for(const auto& J : data) data_to_write.col(idx++) = J;

            group.create_dataset("Element " + std::to_string(I->get_tag()), data_to_write.mem, {static_cast<hsize_t>(data_to_write.n_cols), static_cast<hsize_t>(data_to_write.n_rows)});
        }
    }
#endif
}

FrameRecorder::FrameRecorder(const unsigned T, const OutputType L, const unsigned I)
    : Recorder(T, {}, L, I, true) {
#ifdef SUANPAN_HDF5
    std::ostringstream file_name;
    file_name << 'R' << get_tag() << '-' << to_name(original_type) << ".h5";
    file = hdf5::File(SUANPAN_OUTPUT / file_name.str(), hdf5::File::Mode::Truncate);
#endif
}

void FrameRecorder::save(const shared_ptr<DomainBase>&) {}

void FrameRecorder::print() { suanpan_info("A frame recorder.\n"); }
