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
 * @class Rotation2D
 * @brief A Rotation2D material class.
 * @author tlc
 * @date 10/10/2019
 * @version 0.1.1
 * @file Rotation2D.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef ROTATION2D_H
#define ROTATION2D_H

#include <Material/Material2D/Material2D.h>

class Rotation2D final : public Material2D {
    const unsigned mat_tag;

    unique_ptr<Material> mat_obj;

    mat trans_mat;

public:
    Rotation2D(unsigned, // tag
               unsigned, // mat tag
               double    // Euler angle
    );
    Rotation2D(const Rotation2D&);
    Rotation2D(Rotation2D&&) = delete;
    Rotation2D& operator=(const Rotation2D&) = delete;
    Rotation2D& operator=(Rotation2D&&) = delete;
    ~Rotation2D() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
