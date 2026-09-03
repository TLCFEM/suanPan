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
 * @class NonlocalIsotropicElastic3D
 * @brief The NonlocalIsotropicElastic3D class defines a isotropic elastic material for 3-D
 * problems.
 *
 * The Young's modulus is stored in `elastic_modulus`.
 * The Poisson's ratio is stored in `poissons_ratio`.
 *
 * algorithm verified at 24 April 2019 by tlc
 *
 * @author tlc
 * @date 14/08/2026
 * @version 1.0.0
 * @file NonlocalIsotropicElastic3D.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLOCALISOTROPICELASTIC3D_H
#define NONLOCALISOTROPICELASTIC3D_H

#include <Material/Material3D/Material3D.h>

struct DataNonlocalIsotropicElastic3D {
    const double elastic_modulus;
    const double poissons_ratio;
    const double maximum_energy;
    const double evolution_rate;
    const double reference_length;
    const double diffusion_rate;
};

class NonlocalIsotropicElastic3D final : protected DataNonlocalIsotropicElastic3D, public NonlocalMaterial3D {
public:
    enum class EnergyType {
        TOTAL,
        TENSILE,
        DEV_TENSILE
    };

private:
    static const uvec u_dof, d_dof;

    const EnergyType energy_type;

    [[nodiscard]] std::pair<double, double> compute_scale(const double k) const {
        const auto exp_term = std::exp(diffusion_rate * k);
        const auto s = exp_term / reference_length;

        return {s, diffusion_rate * s};
    }

public:
    NonlocalIsotropicElastic3D(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // poissons ratio
        double,     // maximum energy
        double,     // evolution rate
        double,     // reference length
        double,     // evolution rate
        EnergyType, // energy type
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get(Parameter) const override;

    unique_ptr<Material> unique_copy() override;

    [[nodiscard]] unsigned nonlocal_size() const override;

    int update_trial_status(const vec&) override;

    void print() override;
};

#endif

//! @}
