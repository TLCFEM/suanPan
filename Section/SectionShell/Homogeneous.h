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
 * @class Homogeneous
 * @brief A Homogeneous class.
 * @author tlc
 * @date 18/10/2019
 * @version 0.1.0
 * @file Homogeneous.h
 * @addtogroup SectionShell
 * @{
 */

#ifndef HOMOGENEOUS_H
#define HOMOGENEOUS_H

#include <Section/SectionShell/SectionShell.h>

class Homogeneous final : public SectionShell {
    const unsigned num_ip;
    const double thickness;

    struct IntegrationPoint final {
        const double eccentricity, factor;
        unique_ptr<Material> s_material;
        IntegrationPoint(double, double, unique_ptr<Material>&&);
        IntegrationPoint(const IntegrationPoint&);
        IntegrationPoint(IntegrationPoint&&) noexcept = default;
        IntegrationPoint& operator=(const IntegrationPoint&) = delete;
        IntegrationPoint& operator=(IntegrationPoint&&) noexcept = delete;
        ~IntegrationPoint() = default;
    };

    vector<IntegrationPoint> int_pt;

public:
    Homogeneous(unsigned,    // unique tag
                unsigned,    // material tag
                double,      // thickness
                unsigned = 5 // number of IPs
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<SectionShell> get_copy() override;

    int update_trial_status(const vec&, const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
