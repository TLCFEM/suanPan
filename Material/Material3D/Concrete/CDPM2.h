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
 * @class CDPM2
 * @brief The CDPM2 class.
 *
 * A 3D concrete material model that supports stiffness degradation.
 *
 * @author tlc
 * @date 04/08/2021
 * @version 1.0.0
 * @file CDPM2.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CDPM2_H
#define CDPM2_H

#include <Material/Material3D/Material3D.h>

struct DataCDPM2 {
    const double elastic_modulus = 3E4;
    const double poissons_ratio = .3;
    const double ft = 3.;
    const double fc = 10.;
    const double qh0 = .3;
    const double hp = .01;
    const double df = .85;
    const double ah = .08;
    const double bh = .003;
    const double ch = 2.;
    const double dh = 1E-6;
    const double as = 5.;
    const double eft = 2E-4;
    const double efc = 1E-4;

    const double e = 1.;
    const double e0 = ft / elastic_modulus;
    const double ftfc = ft / fc;
    const double m0 = 3. * (fc / ft - ftfc) * e / (1. + e);
    const double lndf = log(df + 1.) - log(2. * df - 1.);
    const double sqrtdf = ft * sqrt(2. / (3. + 6. * df * df));
    const double eh = bh - dh;
    const double fh = ch * eh / (ah - bh);
};

class CDPM2 final : protected DataCDPM2, public Material3D {
public:
    enum class DamageType {
        NODAMAGE,
        ISOTROPIC,
        ANISOTROPIC
    };

private:
    static constexpr unsigned max_iteration = 20;
    static const double sqrt_six;
    static const double sqrt_three_two;
    static const mat unit_dev_tensor;

    const double double_shear = elastic_modulus / (1. + poissons_ratio);
    const double bulk = elastic_modulus / (3. - 6. * poissons_ratio);

    const DamageType damage_type = DamageType::ANISOTROPIC;

    void compute_plasticity(double, double, double, podarray<double>&) const;
    int compute_damage(double, double, double, double, double, podarray<double>&);
    int compute_damage_factor(double, double, double, double, double&, podarray<double>&) const;

public:
    CDPM2(unsigned,   // tag
          double,     // elastic_modulus
          double,     // poissons_ratio
          double,     // ft
          double,     // fc
          double,     // qh0
          double,     // hp
          double,     // df
          double,     // ah
          double,     // bh
          double,     // ch
          double,     // dh
          double,     // as
          double,     // eft
          double,     // efc
          DamageType, // damage type
          double      // density
        );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

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
