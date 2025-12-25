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
#ifdef SUANPAN_VTK
#include <vtkVertex.h>
#endif

class Particle : public Element {
public:
    Particle(
        unsigned,                // tag
        unsigned,                // node tag
        std::vector<Node::DOF>&& // dof identifier
    );

    [[nodiscard]] Type type() const final { return Type::DEM; }

    int update_status() override { return SUANPAN_SUCCESS; }

    int commit_status() override { return SUANPAN_SUCCESS; }
    int clear_status() override { return SUANPAN_SUCCESS; }
    int reset_status() override { return SUANPAN_SUCCESS; }
};

template<unsigned DIM> class SphericalParticle : public Particle {
protected:
    const double radius, elastic_modulus, poisson_ratio, damping, mass, inertia;

public:
    SphericalParticle(const unsigned T, const unsigned N, std::vector<Node::DOF>&& D, const double R, const double E, const double V, const double A, const double M, const double I)
        : Particle(T, N, std::move(D))
        , radius(R)
        , elastic_modulus(E)
        , poisson_ratio(V)
        , damping(A)
        , mass(M)
        , inertia(I) {}

    [[nodiscard]] double get(const Parameter P) const override {
        if(Parameter::ELASTIC == P) return elastic_modulus;
        if(Parameter::POISSON == P) return poisson_ratio;
        if(Parameter::RADIUS == P) return radius;
        if(Parameter::MASS == P) return mass;
        if(Parameter::INERTIA == P) return inertia;
        if(Parameter::DAMPING == P) return damping;

        return 0.;
    }

#ifdef SUANPAN_VTK
    [[nodiscard]] vtkSmartPointer<vtkCell> GetCell() const override { return vtkSmartPointer<vtkVertex>::New(); }

    mat GetData(const OutputType P) override {
        if(OutputType::A == P) return get_current_acceleration();
        if(OutputType::V == P) return get_current_velocity();
        if(OutputType::U == P) return get_current_displacement();

        return {};
    }

    mat GetDeformation(double) override { return get_coordinate(DIM).t() + get_current_displacement().head(DIM); }
#endif
};

class InertialSphericalParticle2D final : public SphericalParticle<2u> {
public:
    InertialSphericalParticle2D(
        unsigned, // tag
        unsigned, // node tag
        double,   // radius
        double,   // elastic modulus
        double,   // poisson ratio
        double,   // damping
        double,   // mass
        double    // inertia
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

class SphericalParticle2D final : public SphericalParticle<2u> {
public:
    SphericalParticle2D(
        unsigned, // tag
        unsigned, // node tag
        double,   // radius
        double,   // elastic modulus
        double,   // poisson ratio
        double,   // damping
        double    // mass
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

class InertialSphericalParticle3D final : public SphericalParticle<3u> {
public:
    InertialSphericalParticle3D(
        unsigned, // tag
        unsigned, // node tag
        double,   // radius
        double,   // elastic modulus
        double,   // poisson ratio
        double,   // damping
        double,   // mass
        double    // inertia
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

class SphericalParticle3D final : public SphericalParticle<3u> {
public:
    SphericalParticle3D(
        unsigned, // tag
        unsigned, // node tag
        double,   // radius
        double,   // elastic modulus
        double,   // poisson ratio
        double,   // damping
        double    // mass
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
