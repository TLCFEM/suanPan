/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
#ifdef SUANPAN_GCC
#include <sys/stat.h>
#endif

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

    struct stat buffer{};
    if(mat ext_data; stat(file_name.c_str(), &buffer) == 0 && ext_data.load(file_name, raw_ascii)) {
        if(2 == ext_data.n_cols) {
            time = ext_data.col(0);
            magnitude = ext_data.col(1);
        }
        else if(ext_data.n_cols > 2) suanpan_warning("Tabular() reads more than two columns from the given file, check it.\n");
        else {
            suanpan_error("Tabular() requires two valid columns.\n");
            D->disable_amplitude(get_tag());
        }
    }
    else {
        suanpan_error("cannot load file.\n");
        D->disable_amplitude(get_tag());
    }
}

double Tabular::get_amplitude(const double T) {
    const auto step_time = T - start_time;

    uword IDX = 0;
    while(IDX < time.n_elem && time(IDX) < step_time) ++IDX;

    return IDX == 0 ? 0. : IDX == time.n_elem ? magnitude(magnitude.n_elem - 1) : magnitude(IDX - 1) + (step_time - time(IDX - 1)) * (magnitude(IDX) - magnitude(IDX - 1)) / (time(IDX) - time(IDX - 1));
}

void Tabular::print() { suanpan_info("Tabular.\n"); }
