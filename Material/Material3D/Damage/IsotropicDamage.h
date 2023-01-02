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
 * @class IsotropicDamage
 * @brief A IsotropicDamage material class.
 * @author tlc
 * @date 09/06/2020
 * @version 0.1.0
 * @file IsotropicDamage.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef ISOTROPICDAMAGE_H
#define ISOTROPICDAMAGE_H

#include <Material/Material3D/Material3D.h>

class IsotropicDamage : public Material3D {
    const unsigned mat_tag;

    unique_ptr<Material> mat_ptr;

protected:
    virtual void compute_damage() = 0;

public:
    IsotropicDamage(unsigned, // tag
                    unsigned  // mat tag
        );
    IsotropicDamage(const IsotropicDamage&);
    IsotropicDamage(IsotropicDamage&&) = delete;
    IsotropicDamage& operator=(const IsotropicDamage&) = delete;
    IsotropicDamage& operator=(IsotropicDamage&&) = delete;
    ~IsotropicDamage() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
