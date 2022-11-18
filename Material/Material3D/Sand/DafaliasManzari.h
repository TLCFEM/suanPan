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
 * @class DafaliasManzari
 * @brief The DafaliasManzari class.
 *
 * @author tlc
 * @date 10/07/2021
 * @version 0.1.0
 * @file DafaliasManzari.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef DAFALIASMANZARI_H
#define DAFALIASMANZARI_H

#include <Material/Material3D/Material3D.h>

struct DataDafaliasManzari {
    const double shear_modulus = 125.; // shear modulus
    const double poissons_ratio = .05; // poisson's ratio
    const double ac = 1.25;
    const double lc = .02;
    const double e0 = .9;
    const double xi = .7;
    const double m = .01;
    const double h0 = 7.;
    const double h1 = .1;
    const double ch = .9;
    const double nb = 1.1;
    const double a = -.7;
    const double nd = 3.5;
    const double zm = 4.;
    const double cz = 6E2;
    const double pc = -130.;
    const double gr = .2;
};

class DafaliasManzari final : DataDafaliasManzari, public Material3D {
    static constexpr unsigned max_iteration = 20;
    static const mat66 unit_dev_tensor;

    const double pr = (2. + 2. * poissons_ratio) / (3. - 6. * poissons_ratio);
    const double gi = gr * shear_modulus * abs(pc);

    static constexpr uword sa = 0, si = 0, sj = 1;
    static const span sb, sk, sl, sm;

public:
    DafaliasManzari(unsigned,   // tag
                    double,     // g0
                    double,     // nu
                    double,     // ac
                    double,     // lc
                    double,     // e0
                    double,     // xi
                    double,     // m
                    double,     // h0
                    double,     // h1
                    double,     // ch
                    double,     // nb
                    double,     // a
                    double,     // nd
                    double,     // zm
                    double,     // cz
                    double,     // pc
                    double,     // gr
                    double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
