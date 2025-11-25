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
 * @class Recorder
 * @brief A Recorder class.
 * @author tlc
 * @date 27/07/2017
 * @version 0.1.0
 * @file Recorder.h
 * @{
 */

#ifndef RECORDER_H
#define RECORDER_H

#include <Domain/Tag.h>
#include <Recorder/OutputType.h>

class DomainBase;

class Recorder : public UniqueTag {
    const bool use_hdf5;

    static auto normalise_size(std::vector<std::vector<vec>>&);

protected:
    static std::vector<vec> normalise_size(std::vector<vec>&&);

    const OutputType original_type, variable_type;
    const int component;
    const unsigned interval;
    const uvec reference_tag;

    uvec object_tag;
    std::vector<double> time_pool;                        // recorded data
    std::vector<std::vector<std::vector<vec>>> data_pool; // recorded data

    unsigned counter = 0u;

    virtual const uvec& update_tag(const shared_ptr<DomainBase>&);

    virtual void record_impl(const shared_ptr<DomainBase>&) = 0;

public:
    Recorder(
        unsigned,   // tag
        uvec&&,     // object tags
        OutputType, // recorder type
        unsigned,   // interval
        bool        // if to use hdf5
    );

    virtual void initialize(const shared_ptr<DomainBase>&);

    void insert(double);
    void insert(std::vector<vec>&&, unsigned);

    void record(const shared_ptr<DomainBase>&);

    virtual void clear_status();

    virtual void save();

    void print() override;
};

#endif

//! @}
