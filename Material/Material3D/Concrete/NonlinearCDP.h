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
 * @class NonlinearCDP
 * @brief The NonlinearCDP class.
 *
 * A 3D concrete material model that supports stiffness degradation.
 *
 * algorithm verified at 29 April 2019 by tlc
 *
 * References:
 *     1. A Plastic-Damage Model for Concrete.
 *     https://doi.org/10.1016/0020-7683(89)90050-4
 *     2. Plastic-Damage Model for Cyclic Loading of Concrete Structures.
 *     https://doi.org/10.1061/(ASCE)0733-9399(1998)124:8(892)
 *     3. A Plastic-Damage Concrete Model for Earthquake Analysis of Dams.
 *     https://doi.org/10.1002/(SICI)1096-9845(199809)27:9<937::AID-EQE764>3.0.CO;2-5
 *     4. A Return-Mapping Algorithm for Plastic-Damage Models: 3-D and Plane Stress Formulation.
 *     https://doi.org/10.1002/1097-0207(20010120)50:2<487::AID-NME44>3.0.CO;2-N
 *
 * @author tlc
 * @date 29/04/2019
 * @version 1.0.0
 * @file NonlinearCDP.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARCDP_H
#define NONLINEARCDP_H

#include <Material/Material3D/Material3D.h>

struct DataNonlinearCDP {
    const double elastic_modulus, poissons_ratio;
    const double g_t, g_c;
    const double alpha, alpha_p, s0;
};

class NonlinearCDP : protected DataNonlinearCDP, public Material3D {
    static constexpr unsigned max_iteration = 20;
    static const double root_three_two;
    static constexpr double scale = .999; // to avoid overshoot
    static const mat unit_dev_tensor;

    const double bulk = elastic_modulus / (3. - 6. * poissons_ratio);
    const double double_shear = elastic_modulus / (1. + poissons_ratio);
    const double three_alpha_p_bulk = 3. * alpha_p * bulk;
    const double pfplambda = -3. * alpha * three_alpha_p_bulk - root_three_two * double_shear;
    const double one_minus_alpha = 1. - alpha;

    const vec unit_alpha_p{alpha_p, alpha_p, alpha_p, 0., 0., 0.};

    static double compute_r(const vec&);
    static vec compute_dr(const vec&);

    [[nodiscard]] inline double compute_s(double) const;

    /**
     * \brief compute tension backbone
     *
     * return a podarray that contains six parameters
     * 1. d: the damage index
     * 2. f: the final stress = f=(1-d)\bar{f}
     * 3. \bar{f}: the effective stress = f/(1-d)
     * 4. \md{d}: derivative of d
     * 5. \md{f}: derivative of f
     * 6. \md{\bar{f}}: derivative of \bar{f}
     *
     * \return d f \bar{f} \md{d} \md{f} \md{\bar{f}}
     */
    [[nodiscard]] virtual podarray<double> compute_tension_backbone(double) const = 0;

    /**
     * \brief compute compression backbone
     *
     * return a podarray that contains six parameters
     * 1. d: the damage index
     * 2. f: the final stress = f=(1-d)\bar{f}
     * 3. \bar{f}: the effective stress = f/(1-d)
     * 4. \md{d}: derivative of d
     * 5. \md{f}: derivative of f
     * 6. \md{\bar{f}}: derivative of \bar{f}
     *
     * \return d f \bar{f} \md{d} \md{f} \md{\bar{f}}
     */
    [[nodiscard]] virtual podarray<double> compute_compression_backbone(double) const = 0;

public:
    NonlinearCDP(unsigned, // tag
                 double,   // elastic modulus
                 double,   // poissons ratio
                 double,   // normalized crack energy (+)
                 double,   // normalized crush energy (+)
                 double,   // dilatancy parameter
                 double,   // biaxial compression strength ratio
                 double,   // stiffness recovery
                 double    // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
