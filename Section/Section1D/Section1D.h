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
 * @class Section1D
 * @brief A Section1D class.
 * @author tlc
 * @date 27/10/2017
 * @version 0.1.0
 * @file Section1D.h
 * @addtogroup Section-1D
 * @ingroup Section
 * @{
 */

#ifndef SECTION1D_H
#define SECTION1D_H

#include <Section/Section.h>

using std::vector;

class Section1D : public Section {
protected:
    unique_ptr<Material> s_material;

public:
    explicit Section1D(unsigned, unsigned, double);
    Section1D(const Section1D&);
    Section1D(Section1D&&) = delete;                 // move forbidden
    Section1D& operator=(const Section1D&) = delete; // assign forbidden
    Section1D& operator=(Section1D&&) = delete;      // assign forbidden
    ~Section1D() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
