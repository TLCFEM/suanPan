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
/**
 * @class FrameRecorder
 * @brief A FrameRecorder class.
 *
 * @author tlc
 * @date 30/03/2019
 * @version 0.1.0
 * @file FrameRecorder.h
 * @addtogroup Recorder
 * @{
 */

#ifndef FRAMERECORDER_H
#define FRAMERECORDER_H

#include <Recorder/Recorder.h>

class FrameRecorder final : public Recorder {
#ifdef SUANPAN_HDF5
    hid_t file_id = 0;
#endif

public:
    FrameRecorder(
        unsigned,   // tag
        OutputType, // recorder type
        unsigned    // interval
    );
    FrameRecorder(const FrameRecorder&) = delete;
    FrameRecorder(FrameRecorder&&) = delete;
    FrameRecorder& operator=(const FrameRecorder&) = delete;
    FrameRecorder& operator=(FrameRecorder&&) = delete;
    ~FrameRecorder() override;

    void record(const shared_ptr<DomainBase>&) override;

    void save() override;

    void print() override;
};

#endif

//! @}
