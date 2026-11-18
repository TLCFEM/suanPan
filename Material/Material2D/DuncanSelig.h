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
 * @class DuncanSelig
 * @brief A DuncanSelig material class.
 *
 * For plane strain soil.
 *
 * @author tlc
 * @date 01/04/2024
 * @version 1.0.0
 * @file DuncanSelig.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef DUNCANSELIG_H
#define DUNCANSELIG_H

#include <Material/Material2D/Material2D.h>

class DuncanSelig final : public Material2D {
    static constexpr unsigned max_iteration = 30u;
    static constexpr double min_ratio = .1;
    static constexpr double min_degradation = .1;

    static double dev(const vec&);
    static mat compute_stiffness(double, double);

    double p_atm = 14.7;
    double ref_elastic = 400. * p_atm, n = .6;
    double ref_bulk = 300. * p_atm, m = .2;
    double ini_phi = .7, ten_fold_phi_diff = .1, r_f = .7, cohesion = .5;

    using ds_tuple = std::tuple<double, double, rowvec3, rowvec3>;

    [[nodiscard]] std::tuple<double, double> compute_modulus(double, double, double) const;
    [[nodiscard]] std::tuple<double, double> compute_elastic(double) const;
    [[nodiscard]] std::tuple<double, double> compute_bulk(double) const;
    [[nodiscard]] ds_tuple compute_principal_stress(const vec&) const;
    [[nodiscard]] ds_tuple compute_elastic_moduli(const vec&) const;
    [[nodiscard]] ds_tuple compute_plastic_moduli(const vec&) const;

    template<bool plastic> [[nodiscard]] auto compute_moduli(const vec& t_stress) const {
        if constexpr(plastic) return compute_plastic_moduli(t_stress);
        return compute_elastic_moduli(t_stress);
    }

    template<bool plastic, bool with_jacobian> [[nodiscard]] auto local_assemble(const vec& t_stress) {
        const auto [elastic, bulk, deds, dkds] = compute_moduli<plastic>(t_stress);

        const auto factor_a = 3. * bulk * (3. * bulk + elastic);
        const auto factor_b = 3. * bulk * (3. * bulk - elastic);

        mat33 right(fill::zeros);
        right(0, 0) = right(1, 1) = factor_a;
        right(0, 1) = right(1, 0) = factor_b;
        right(2, 2) = 3. * bulk * elastic;

        const vec3 stress_diff = t_stress - current_stress;
        const vec3 residual = (9. * bulk - elastic) * stress_diff - right * incre_strain;

        mat33 jacobian;

        if constexpr(with_jacobian) {
            const rowvec3 factor_c = 3. * (elastic * dkds + bulk * deds);
            const rowvec3 dfads = 18. * bulk * dkds + factor_c;
            const rowvec3 dfbds = 18. * bulk * dkds - factor_c;

            jacobian = (9. * bulk - elastic) * eye(3, 3) + stress_diff * (9. * dkds - deds);
            jacobian.row(0) -= incre_strain(0) * dfads + incre_strain(1) * dfbds;
            jacobian.row(1) -= incre_strain(0) * dfbds + incre_strain(1) * dfads;
            jacobian.row(2) -= incre_strain(2) * factor_c;
        }

        return std::tuple{residual, right, jacobian};
    };

    template<bool plastic> int local_update() {
        auto ref_error = 0.;
        vec3 incre;

        auto counter = 0u;
        while(true) {
            if(max_iteration < ++counter) {
                suanpan_debug("Local iteration cannot converge within {} iterations.\n", max_iteration);
                return SUANPAN_FAIL;
            }

            const auto [residual, right, jacobian] = local_assemble<plastic, true>(trial_stress);

            if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

            const auto error = suanpan::inf_norm(incre);
            if(1u == counter) ref_error = error;
            suanpan_debug("Local iteration error: {:.5E}.\n", error);
            const auto residual_norm = suanpan::inf_norm(residual);
            if(error < tolerance * ref_error || ((error < tolerance || residual_norm < tolerance) && counter > 5u)) {
                if(!solve(trial_stiffness, jacobian, right)) return SUANPAN_FAIL;

                return SUANPAN_SUCCESS;
            }

            // backtracking line search on the residual norm
            auto step = 1.;
            vec candidate = trial_stress - incre;
            while(step > 1E-3 && suanpan::inf_norm(std::get<0>(local_assemble<plastic, false>(candidate))) >= residual_norm) candidate = trial_stress - (step *= .5) * incre;

            trial_stress = candidate;
        }
    }

public:
    DuncanSelig(
        unsigned,   // tag
        const vec&, // parameters
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> unique_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
