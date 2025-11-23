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
 * @class Particle
 * @brief The Particle class represents the abstract particle used in DEM.
 *
 * @author tlc
 * @date 23/11/2025
 * @version 0.1.0
 * @file Particle.h
 * @addtogroup Particle
 * @ingroup Element
 * @{
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <Element/Element.h>

enum class ParticleParameter {
    ELASTICMODULUS,
    POISSONSRATIO,
    RADIUS,
    MASS,
    INERTIA
};

class Particle : public Element {
public:
    Particle(
        unsigned,          // tag
        unsigned,          // node tag
        std::vector<DOF>&& // dof identifier
    );

    int update_status() override { return SUANPAN_SUCCESS; }

    int commit_status() override { return SUANPAN_SUCCESS; }
    int clear_status() override { return SUANPAN_SUCCESS; }
    int reset_status() override { return SUANPAN_SUCCESS; }
};

class SphericalParticle : public Particle {
    const double radius, elastic_modulus, poisson_ratio, mass, inertia;

public:
    SphericalParticle(
        unsigned,           // tag
        unsigned,           // node tag
        std::vector<DOF>&&, // dof identifier
        double,             // radius
        double,             // elastic modulus
        double,             // poisson ratio
        double,             // mass
        double              // inertia
    );

    [[nodiscard]] double get_parameter(ParticleParameter) const;
};

class InertialSphericalParticle2D : public SphericalParticle {
public:
    InertialSphericalParticle2D(
        unsigned, // tag
        unsigned, // node tag
        double,   // radius
        double,   // elastic modulus
        double,   // poisson ratio
        double,   // mass
        double    // inertia
    );
};

class SphericalParticle2D : public SphericalParticle {
public:
    SphericalParticle2D(
        unsigned, // tag
        unsigned, // node tag
        double,   // radius
        double,   // elastic modulus
        double,   // poisson ratio
        double    // mass
    );
};

#endif

//! @}
