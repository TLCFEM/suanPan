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
 * @class SectionNM
 * @brief A SectionNM class.
 * @author tlc
 * @date 26/11/2021
 * @version 0.1.0
 * @file SectionNM.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef SECTIONNM_H
#define SECTIONNM_H

#include <Section/Section.h>

struct DataSectionNM {
    vec initial_history;
    vec current_history;
    vec trial_history;
};

class SectionNM : protected DataSectionNM, public Section {
protected:
    static constexpr double tolerance = 1E-12;

public:
    using Section::Section;

    void initialize_history(unsigned);

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
