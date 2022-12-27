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

#include "AmplitudeRecorder.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/Amplitude/Amplitude.h>

void AmplitudeRecorder::initialize(const shared_ptr<DomainBase>& D) {
    for(const auto I : get_object_tag())
        if(!D->find<Amplitude>(I)) {
            D->disable_recorder(get_tag());
            return;
        }
}

void AmplitudeRecorder::record(const shared_ptr<DomainBase>& D) {
    if(!if_perform_record()) return;

    const sp_d auto current_time = D->get_factory()->get_current_time();
    auto& obj_tag = get_object_tag();

    for(unsigned I = 0; I < obj_tag.n_elem; ++I) insert({{D->get<Amplitude>(obj_tag(I))->get_amplitude(current_time)}}, I);

    if(if_record_time()) insert(current_time);
}

void AmplitudeRecorder::print() { suanpan_info("A recorder to record amplitudes.\n"); }
