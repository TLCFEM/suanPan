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

#include "NZStrongMotion.h"

#include <Domain/DomainBase.h>

NZStrongMotion::NZStrongMotion(const unsigned T, const char* P)
    : Amplitude(T)
    , file_name(P) {}

unique_ptr<Amplitude> NZStrongMotion::get_copy() { return std::make_unique<NZStrongMotion>(*this); }

void NZStrongMotion::initialize(const shared_ptr<DomainBase>& D) {
    std::error_code code;
    if(Col<int> data; !fs::exists(file_name, code) || !data.load(file_name, auto_detect)) {
        suanpan_error("Cannot load \"{}\".\n", file_name);
        D->disable_amplitude(get_tag());
    }
    else magnitude = conv_to<vec>::from(data);
}

double NZStrongMotion::get_amplitude(const double T) {
    const auto step_time = T - start_time;

    const auto IDX = static_cast<uword>(step_time * 1E3 / magnitude(0));

    if(IDX >= magnitude.n_elem - 3) return 0.;

    const auto diff_magnitude = magnitude(IDX + 3) - magnitude(IDX + 2);

    return (magnitude(IDX + 2) + (step_time - static_cast<double>(IDX) * magnitude(0) / 1E3) * diff_magnitude / magnitude(0)) / magnitude(1);
}

void NZStrongMotion::print() {
    suanpan_info("NZStrongMotion {} with total duration of {:.2f} seconds and raw PGA of {:.4f} mm/s/s.\n", file_name, static_cast<double>(magnitude.n_elem - 2llu) * magnitude(0) / 1E3, magnitude(1) / 1E3);
}
