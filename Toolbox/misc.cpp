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

#include "misc.h"

void save_result(const mat& result) {
#ifdef SUANPAN_HDF5
    if(!result.save("RESULT.h5", hdf5_binary_trans))
        suanpan_error("Fail to save to file.\n");
#else
    if(!result.save("RESULT.txt", raw_ascii))
        suanpan_error("Fail to save to file.\n");
#endif
}

void save_gnuplot() {
    if(std::ofstream gnuplot("RESULT.plt"); gnuplot.is_open()) {
        gnuplot << "reset\n";
        gnuplot << "set term tikz size 14cm,10cm\n";
        gnuplot << "set output \"RESULT.tex\"\n";
        gnuplot << "unset key\n";
        gnuplot << "set xrange [*:*]\n";
        gnuplot << "set yrange [*:*]\n";
        gnuplot << "set xlabel \"input\"\n";
        gnuplot << "set ylabel \"output\"\n";
        gnuplot << "set grid\n";
        gnuplot << "plot \"RESULT.txt\" u 1:2 w l lw 2\n";
        gnuplot << "set output\n";
    }
}
