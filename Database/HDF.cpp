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

#ifdef SUANPAN_HDF5

#include "HDF.h"
#include <hdf5/hdf5.h>

//#include <hdf5_hl.h>

HDF::HDF(const char* N)
    : file_id(H5Fcreate(N, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) {}

HDF::~HDF() { H5Fclose(file_id); }

int HDF::save() { return 0; }

#endif
