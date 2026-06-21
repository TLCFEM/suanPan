/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

Particle::Particle(const unsigned T, const unsigned N, std::vector<Node::DOF>&& D)
    : Element(T, 1, static_cast<unsigned>(D.size()), uvec{N}, std::move(D)) {}

InertialSphericalParticle2D::InertialSphericalParticle2D(const unsigned T, const unsigned N, const double R, const double E, const double V, const double A, const double M, const double I)
    : SphericalParticle(T, N, {Node::DOF::U1, Node::DOF::U2, Node::DOF::UR3}, R, E, V, A, M, I) {}

SP_STATUS InertialSphericalParticle2D::initialize(const shared_ptr<DomainBase>&) {
    initial_mass = diagmat(vec{mass, mass, inertia});

    ConstantMass(this);

    return SP_STATUS::SUCCESS;
}

SphericalParticle2D::SphericalParticle2D(const unsigned T, const unsigned N, const double R, const double E, const double V, const double A, const double M)
    : SphericalParticle(T, N, {Node::DOF::U1, Node::DOF::U2}, R, E, V, A, M, 0.) {}

SP_STATUS SphericalParticle2D::initialize(const shared_ptr<DomainBase>&) {
    initial_mass = diagmat(vec{mass, mass});

    ConstantMass(this);

    return SP_STATUS::SUCCESS;
}

InertialSphericalParticle3D::InertialSphericalParticle3D(const unsigned T, const unsigned N, const double R, const double E, const double V, const double A, const double M, const double I)
    : SphericalParticle(T, N, {Node::DOF::U1, Node::DOF::U2, Node::DOF::U3, Node::DOF::UR1, Node::DOF::UR2, Node::DOF::UR3}, R, E, V, A, M, I) {}

SP_STATUS InertialSphericalParticle3D::initialize(const shared_ptr<DomainBase>&) {
    initial_mass = diagmat(vec{mass, mass, mass, inertia, inertia, inertia});

    ConstantMass(this);

    return SP_STATUS::SUCCESS;
}

SphericalParticle3D::SphericalParticle3D(const unsigned T, const unsigned N, const double R, const double E, const double V, const double A, const double M)
    : SphericalParticle(T, N, {Node::DOF::U1, Node::DOF::U2, Node::DOF::U3}, R, E, V, A, M, 0.) {}

SP_STATUS SphericalParticle3D::initialize(const shared_ptr<DomainBase>&) {
    initial_mass = diagmat(vec{mass, mass, mass});

    ConstantMass(this);

    return SP_STATUS::SUCCESS;
}
