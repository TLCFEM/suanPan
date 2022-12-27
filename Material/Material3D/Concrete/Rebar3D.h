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
 * @class Rebar3D
 * @brief A Rebar3D material class.
 * @author tlc
 * @date 11/12/2017
 * @version 0.1.1
 * @file Rebar3D.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef REBAR3D_H
#define REBAR3D_H

#include <Material/Material3D/Material3D.h>

class Rebar3D final : public Material3D {
    const unsigned tag_x, tag_y, tag_z;

    const double ratio_x, ratio_y, ratio_z;

    const mat trans_mat;

    unique_ptr<Material> rebar_x, rebar_y, rebar_z;

public:
    Rebar3D(unsigned,   // tag
            unsigned,   // material tag along x axis
            unsigned,   // material tag along y axis
            unsigned,   // material tag along z axis
            double,     // reinforcement ratio along x axis
            double,     // reinforcement ratio along y axis
            double,     // reinforcement ratio along z axis
            double = 0. // inclination
    );
    Rebar3D(const Rebar3D&);
    Rebar3D(Rebar3D&&) = delete;
    Rebar3D& operator=(const Rebar3D&) = delete;
    Rebar3D& operator=(Rebar3D&&) = delete;
    ~Rebar3D() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;
};

#endif

//! @}
