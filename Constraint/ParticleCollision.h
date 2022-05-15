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
 * @class ParticleCollision
 * @brief A ParticleCollision class.
 *
 * @author tlc
 * @date 25/08/2020
 * @version 0.1.0
 * @file ParticleCollision.h
 * @addtogroup Constraint
 * @{
 */

#ifndef PARTICLECOLLISION_H
#define PARTICLECOLLISION_H

#include <Constraint/Constraint.h>
#include <Domain/Node.h>

class ParticleCollision : public Constraint {
    const unsigned num_dof;

    /**
     * \brief The force between two particles is assumed to be a function of distance.
     * \return force
     */
    [[nodiscard]] virtual double compute_f(double) const = 0;
    /**
     * \brief Compute the derivative of force w.r.t. distance between two particles.
     * \return the derivative of force w.r.t. distance
     */
    [[nodiscard]] virtual double compute_df(double) const = 0;

protected:
    [[nodiscard]] vec get_position(const shared_ptr<Node>&) const;

    void apply_contact(const shared_ptr<DomainBase>&, const shared_ptr<Node>&, const shared_ptr<Node>&, bool);

    virtual int process_meta(const shared_ptr<DomainBase>&, bool) = 0;

public:
    ParticleCollision(unsigned, unsigned, unsigned);

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;
    int process_resistance(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
