/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
 * @class Load
 * @brief A Load class.
 *
 * The Load class is in charge of returning load level according to given time
 * increment.
 *
 * @author tlc
 * @date 01/10/2017
 * @version 0.2.0
 * @file Load.h
 * @addtogroup Load
 * @{
 */

#ifndef LOAD_H
#define LOAD_H

#include <Domain/ConditionalModifier.h>

class Load : public ConditionalModifier {
protected:
    static double multiplier;

    const bool mpdc_flag = false;

    const double pattern;

    vec trial_load;
    vec trial_settlement;
    sp_vec reference_load;

    friend void set_load_multiplier(double);

public:
    Load(
        unsigned, // tag
        unsigned, // amplitude tag
        uvec&&,   // node tag
        uvec&&,   // dof tag
        double    // nominal magnitude
    );

    void enable_displacement_control() const;
    [[nodiscard]] bool if_displacement_control() const;

    [[nodiscard]] const vec& get_trial_load() const;
    [[nodiscard]] const vec& get_trial_settlement() const;
    [[nodiscard]] const sp_vec& get_reference_load() const;
};

void set_load_multiplier(double);

class GroupLoad {
    const uvec groups;

protected:
    [[nodiscard]] uvec update_object_tag(const shared_ptr<DomainBase>&) const;

public:
    explicit GroupLoad(uvec&&);
};

#endif

//! @}
