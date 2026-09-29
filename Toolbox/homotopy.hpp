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

#ifndef HOMOTOPY_HPP
#define HOMOTOPY_HPP

#include <suanPan.h>

template<typename JT, typename FT, typename ST> requires is_arma_mat<double, JT> && is_arma_mat<double, FT> && std::invocable<ST, const FT&> && std::same_as<std::invoke_result_t<ST, const FT&>, std::pair<FT, JT>> int homotopy_solve(FT& x, ST&& system, double incre_t, const double tolerance, const unsigned max_evaluation) {
    const auto initial_f = system(x).first;

    constexpr auto min_incre{1e-8};
    constexpr auto max_iteration{20u};

    auto counter{0u};
    auto current_t{0.};
    auto current_x = x;

    while(current_t < 1.) {
        const auto trial_t = std::min(1., current_t + incre_t);
        const FT target_f = (1. - trial_t) * initial_f;
        auto trial_x = current_x;

        FT residual;
        JT jacobian;

        const auto bounding_eval = [&](const FT& in_x) {
            if(++counter >= max_evaluation) return false;
            std::tie(residual, jacobian) = system(in_x);
            residual -= target_f;
            return true;
        };

        auto ref_error{1.};
        auto converged{false};
        auto iteration_used{0u};

        for(auto round{0u}; round < max_iteration; ++round) {
            if(!bounding_eval(trial_x)) return SUANPAN_FAIL;

            FT incre_x;
            if(!solve(incre_x, jacobian, residual, solve_opts::equilibrate)) break;

            const auto error = suanpan::inf_norm(incre_x);
            if(0u == round) ref_error = error;
            suanpan_debug("Homotopy iteration error: {:.5E} at progress {:.3f}.\n", error, trial_t);

            if(error < tolerance * ref_error || ((error < tolerance || suanpan::inf_norm(residual) < tolerance) && round > 5u)) {
                converged = true;
                iteration_used = round;
                break;
            }

            trial_x -= incre_x;
        }

        if(converged) {
            current_t = trial_t;
            current_x = trial_x;

            if(iteration_used <= 3) incre_t = std::min(1. - current_t, incre_t * 1.5);
        }
        else if((incre_t *= .5) < min_incre) return SUANPAN_FAIL;
    }

    x = current_x;

    return SUANPAN_SUCCESS
}

#endif

//! @}
