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
 * @class EigenRecorder
 * @brief A EigenRecorder class.
 *
 * @author tlc
 * @date 06/07/2018
 * @version 0.1.0
 * @file EigenRecorder.h
 * @addtogroup Recorder
 * @{
 */

#ifndef EIGENRECORDER_H
#define EIGENRECORDER_H

#include <Recorder/Recorder.h>

class EigenRecorder final : public Recorder {
    vec eigen_value;
    std::vector<std::map<unsigned, vec>> eigen_pool;

public:
    explicit EigenRecorder(unsigned = 0, bool = true);

    void record(const shared_ptr<DomainBase>&) override;

    void save() override;

    void print() override;
};

#endif

//! @}
