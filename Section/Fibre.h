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
 * @class Fibre
 * @brief A Fibre class.
 * @author tlc
 * @date 14/09/2023
 * @version 0.1.0
 * @file Fibre.h
 * @ingroup Section
 * @{
 */

#ifndef FIBRE_H
#define FIBRE_H

#include <Section/Section.h>
#include <Toolbox/ResourceHolder.h>

class Fibre : public Section {
    const uvec fibre_tag;

protected:
    std::vector<ResourceHolder<Section>> fibre;

public:
    Fibre(unsigned, uvec&&, SectionType);
    Fibre(const Fibre&) = default;
    Fibre(Fibre&&) noexcept = delete;
    Fibre& operator=(const Fibre&) = delete;
    Fibre& operator=(Fibre&&) noexcept = delete;
    ~Fibre() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
