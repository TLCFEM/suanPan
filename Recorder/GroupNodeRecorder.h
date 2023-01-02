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
/**
 * @class GroupNodeRecorder
 * @brief A GroupNodeRecorder class.
 *
 * @author tlc
 * @date 24/08/2017
 * @version 0.1.0
 * @file GroupNodeRecorder.h
 * @addtogroup Recorder
 * @{
 */

#ifndef GROUPNODERECORDER_H
#define GROUPNODERECORDER_H

#include "NodeRecorder.h"

class GroupNodeRecorder final : public NodeRecorder {
    const uvec groups;

    void update_tag(const shared_ptr<DomainBase>&);

public:
    GroupNodeRecorder(unsigned,   // tag
                      uvec&&,     // object tags
                      OutputType, // recorder type
                      unsigned,   // interval
                      bool,       // if to record time
                      bool        // if to use hdf5
        );

    void initialize(const shared_ptr<DomainBase>&) override;

    void print() override;
};

#endif

//! @}
