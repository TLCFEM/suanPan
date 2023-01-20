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
 * @class Fibre1D
 * @brief A Fibre1D class.
 * @author tlc
 * @date 13/10/2017
 * @version 0.1.0
 * @file Fibre1D.h
 * @addtogroup Section-1D
 * @ingroup Section
 * @{
 */

#ifndef FIBRE1D_H
#define FIBRE1D_H

#include <Section/Section1D/Section1D.h>

class Fibre1D final : public Section1D {
    uvec fibre_tag;

    vector<unique_ptr<Section>> fibre;

public:
    Fibre1D(unsigned, uvec&&);
    Fibre1D(const Fibre1D&);
    Fibre1D(Fibre1D&&) = delete;
    Fibre1D& operator=(const Fibre1D&) = delete;
    Fibre1D& operator=(Fibre1D&&) = delete;
    ~Fibre1D() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
