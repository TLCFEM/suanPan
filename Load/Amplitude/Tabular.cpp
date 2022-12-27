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

#include "Tabular.h"
#include <Domain/DomainBase.h>

Tabular::Tabular(const unsigned T, vec&& TI, vec&& M, const unsigned ST)
    : Amplitude(T, ST)
    , time(std::forward<vec>(TI))
    , magnitude(std::forward<vec>(M)) {}

Tabular::Tabular(const unsigned T, string&& P, const unsigned ST)
    : Amplitude(T, ST)
    , file_name(std::forward<string>(P)) {}

void Tabular::initialize(const shared_ptr<DomainBase>& D) {
    if(file_name.empty()) {
        if(time.n_elem != magnitude.n_elem) D->disable_amplitude(get_tag());
        return;
    }

    std::error_code code;
    if(mat ext_data; !fs::exists(file_name, code) || !ext_data.load(file_name, raw_ascii)) {
        suanpan_error("Tabular() cannot load file.\n");
        D->disable_amplitude(get_tag());
    }
    else if(ext_data.n_cols >= 2llu) {
        if(ext_data.n_cols > 2llu) suanpan_warning("Tabular() reads more than two columns from the given file, please ensure the correct file is used.\n");
        time = ext_data.col(0);
        magnitude = ext_data.col(1);
    }
    else {
        suanpan_error("Tabular() requires two valid columns.\n");
        D->disable_amplitude(get_tag());
    }
}

double Tabular::get_amplitude(const double T) {
    const auto step_time = T - start_time;

    if(step_time <= time.front()) return magnitude.front();

    if(step_time >= time.back()) return magnitude.back();

    vec result(1);

    interp1(time, magnitude, vec{step_time}, result, "*linear");

    return result(0);
}

void Tabular::print() { suanpan_info("Tabular.\n"); }
