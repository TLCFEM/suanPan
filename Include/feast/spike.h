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

#ifndef SPIKE_H
#define SPIKE_H

#ifdef __cplusplus
extern "C" {
#endif
void spikeinit_(int*, int*, int*);
void dspike_tune_(int*);
void dspike_gbsv_(int*, int*, int*, int*, int*, double*, int*, double*, int*, int*);
void dspike_gbtrf_(int*, int*, int*, int*, double*, int*, double*, int*);
void dspike_gbtrs_(int*, const char*, int*, int*, int*, int*, double*, int*, double*, double*, int*);
void sspike_tune_(int*);
void sspike_gbsv_(int*, int*, int*, int*, int*, float*, int*, float*, int*, int*);
void sspike_gbtrf_(int*, int*, int*, int*, float*, int*, float*, int*);
void sspike_gbtrs_(int*, const char*, int*, int*, int*, int*, float*, int*, float*, float*, int*);
#ifdef __cplusplus
}
#endif

#endif // SPIKE_H
