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

class Recorder : public Tag {
    uvec object_tag;
    OutputType variable_type;
    std::vector<double> time_pool;                        // recorded data
    std::vector<std::vector<std::vector<vec>>> data_pool; // recorded data

    const bool record_time;
    const bool use_hdf5;

protected:
    const unsigned interval;
    unsigned counter = 0;

    bool if_perform_record();

public:
    Recorder(
        unsigned,   // tag
        uvec&&,     // object tags
        OutputType, // recorder type
        unsigned,   // interval
        bool,       // if to record time
        bool        // if to use hdf5
    );
    Recorder(const Recorder&) = delete;
    Recorder(Recorder&&) = delete;                 // move forbidden
    Recorder& operator=(const Recorder&) = delete; // assign forbidden
    Recorder& operator=(Recorder&&) = delete;      // assign forbidden
    ~Recorder() override = default;

    virtual void initialize(const shared_ptr<DomainBase>&);

    void set_object_tag(uvec&&);
    [[nodiscard]] const uvec& get_object_tag() const;

    void set_variable_type(OutputType);
    [[nodiscard]] const OutputType& get_variable_type() const;

    [[nodiscard]] bool if_hdf5() const;
    [[nodiscard]] bool if_record_time() const;

    void insert(double);
    void insert(const std::vector<vec>&, unsigned);

    [[nodiscard]] const std::vector<std::vector<std::vector<vec>>>& get_data_pool() const;
    [[nodiscard]] const std::vector<double>& get_time_pool() const;

    virtual void record(const shared_ptr<DomainBase>&) = 0;

    void clear_status();

    virtual void save();

    void print() override;
};

#endif

//! @}
