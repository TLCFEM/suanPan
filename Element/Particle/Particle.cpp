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

#include "Particle.h"

#include <Domain/DOF.h>

Particle::Particle(const unsigned T, const unsigned N, std::vector<DOF>&& D)
    : Element(T, 1, D.size(), uvec{N}, std::move(D)) {}

SphericalParticle::SphericalParticle(const unsigned T, const unsigned N, std::vector<DOF>&& D, const double R, const double E, const double V, const double M, const double I)
    : Particle(T, N, std::move(D))
    , radius(R)
    , elastic_modulus(E)
    , poisson_ratio(V)
    , mass(M)
    , inertia(I) {}

double SphericalParticle::get_parameter(const ParticleParameter P) const {
    if(ParticleParameter::ELASTICMODULUS == P) return elastic_modulus;
    if(ParticleParameter::POISSONSRATIO == P) return poisson_ratio;
    if(ParticleParameter::RADIUS == P) return radius;
    if(ParticleParameter::MASS == P) return mass;
    if(ParticleParameter::INERTIA == P) return inertia;

    return 0.;
}

InertialSphericalParticle2D::InertialSphericalParticle2D(const unsigned T, const unsigned N, const double R, const double E, const double V, const double M, const double I)
    : SphericalParticle(T, N, {DOF::U1, DOF::U2, DOF::UR3}, R, E, V, M, I) {}

SphericalParticle2D::SphericalParticle2D(const unsigned T, const unsigned N, const double R, const double E, const double V, const double M)
    : SphericalParticle(T, N, {DOF::U1, DOF::U2}, R, E, V, M, 0.) {}
