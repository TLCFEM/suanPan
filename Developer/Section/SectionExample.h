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
 * @class SectionExample
 * @brief A SectionExample class.
 * @author tlc
 * @date 15/09/2017
 * @version 0.1.0
 * @file SectionExample.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef SECTIONEXAMPLE_H
#define SECTIONEXAMPLE_H

#include <Section/Section.h>

class SectionExample final : public Section {
    const double edge_length, moment_inertia, elastic_modulus;

public:
    explicit SectionExample(unsigned, // tag
                            double,   // edge length
                            double    // elastic modulus
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

SUANPAN_EXPORT void new_sectionexample(unique_ptr<Section>&, istringstream&);

#endif

//! @}
