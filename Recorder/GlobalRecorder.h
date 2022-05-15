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
/**
 * @class GlobalRecorder
 * @brief A GlobalRecorder class.
 *
 * @author tlc
 * @date 19/07/2020
 * @version 0.1.0
 * @file GlobalRecorder.h
 * @addtogroup Recorder
 * @{
 */

#ifndef GLOBALRECORDER_H
#define GLOBALRECORDER_H

#include <Recorder/Recorder.h>

class GlobalRecorder : public Recorder {
public:
    GlobalRecorder(unsigned,   // tag
                   OutputType, // recorder type
                   unsigned,   // interval
                   bool,       // if to record time
                   bool        // if to use hdf5
    );

    void record(const shared_ptr<DomainBase>&) override;

    void print() override;
};

#endif

//! @}
