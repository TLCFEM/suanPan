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
 * @class VisualisationRecorder
 * @brief A VisualisationRecorder class.
 *
 * @author tlc
 * @date 22/05/2020
 * @version 0.1.0
 * @file VisualisationRecorder.h
 * @addtogroup Recorder
 * @{
 */

#ifndef VISUALISATIONRECORDER_H
#define VISUALISATIONRECORDER_H

#include <Recorder/Recorder.h>

#ifdef SUANPAN_VTK
#include <Element/Visualisation/vtkParser.h>
#endif

class VisualisationRecorder final : public Recorder {
    unsigned total_counter = 0u;
    unsigned width = 6u;

#ifdef SUANPAN_VTK
    vtkInfo config;

    void (*function_handler)(const shared_ptr<DomainBase>&, vtkInfo) = nullptr;
#endif

public:
    VisualisationRecorder(
        unsigned,   // tag
        OutputType, // recorder type
        unsigned,   // interval
        unsigned,   // output width
        double = 1. // scale
    );

    void record(const shared_ptr<DomainBase>&) override;

    void save() override;

    void print() override;
};

#endif

//! @}
