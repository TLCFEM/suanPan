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
 * @class Interaction
 * @author tlc
 * @date 01/12/2025
 * @version 0.1.0
 * @file Interaction.h
 * @addtogroup Interaction
 * @{
 */

#ifndef INTERACTION_H
#define INTERACTION_H

#include <Domain/Tag.h>

class DomainBase;
class Element;
template<sp_d T> class Factory;

class InteractionPair {
    const shared_ptr<Element> object_i, object_j;

    bool inertial = false;
    unsigned dimension{};

public:
    double effective_mass{};
    double effective_radius{};
    double effective_modulus{};
    double effective_damping{};

    InteractionPair(const shared_ptr<Element>&, const shared_ptr<Element>&);

    void set_dimension(const unsigned dim) { dimension = dim; }
    void set_inertial(const bool flag) { inertial = flag; }

    [[nodiscard]] vec position_i() const;
    [[nodiscard]] vec position_j() const;

    [[nodiscard]] const uvec& dof_i() const;
    [[nodiscard]] const uvec& dof_j() const;

    [[nodiscard]] double initial_gap() const;

    [[nodiscard]] vec relative_velocity() const;
};

class Interaction : public CopyableTag {
protected:
    shared_ptr<Factory<double>> factory;

public:
    using CopyableTag::CopyableTag;

    void initialize(const shared_ptr<DomainBase>&);

    virtual void apply(const shared_ptr<InteractionPair>&) const = 0;
};

class Hertzian final : public Interaction {
    static constexpr double two_third = 2. / 3.;

public:
    using Interaction::Interaction;

    void apply(const shared_ptr<InteractionPair>&) const override;
};

class HertzianDamped final : public Interaction {
    static constexpr double four_third = 4. / 3.;

public:
    using Interaction::Interaction;

    void apply(const shared_ptr<InteractionPair>&) const override;
};

#endif

//! @}
