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
 * @class LBFGS
 * @brief The LBFGS class defines a solver using LBFGS iteration method.
 *
 * The LBFGS method is a rank two quasi--Newton method which has a super linear
 * convergence rate. The LBFGS class supports both conventional BFGS and LBFGS
 * method which uses limited history information.
 *
 * \f{gather}{
 * K_{n+1}^{-1}=\left(I-\dfrac{\Delta{}UR^T}{R^T\Delta{}U}\right)K_n^{-1}\left(I-\dfrac{R\Delta{}U^T}{R^T\Delta{}U}\right)+\dfrac{\Delta{}U\Delta{}U^T}{R^T\Delta{}U}.
 * \f}
 *
 * The \f$I\f$ is identity matrix. The \f$\Delta{}U\f$ is current displacement
 * increment. The \f$R\f$ is current residual. For brevity, in both terms, the
 * subscript \f$n\f$ representing current step is dropped.
 *
 * @author tlc
 * @date 13/07/2022
 * @version 0.1.0
 * @file LBFGS.hpp
 * @addtogroup Utility
 * @{
 */

#ifndef LBFGS_HPP
#define LBFGS_HPP

#include <suanPan.h>
#include <deque>

using std::deque;
using std::vector;

template<typename T> concept Differentiable = requires(T t, const vec& x)
{
    t.evaluate_residual(x);
    t.evaluate_jacobian(x);
};

class LBFGS final {
    deque<vec> hist_ninja, hist_residual;
    deque<double> hist_factor;
    vector<double> alpha;

    const unsigned max_hist;
    const unsigned max_iteration;
    const double abs_tol;
    const double rel_tol;

public:
    explicit LBFGS(const unsigned MH = 10, const unsigned MI = 500, const double AT = 1E-12, const double RT = 1E-12)
        : max_hist(MH)
        , max_iteration(MI)
        , abs_tol(AT)
        , rel_tol(RT) {}

    template<Differentiable F> int optimize(F& func, vec& x) {
        // clear container
        hist_ninja.clear();
        hist_residual.clear();
        hist_factor.clear();

        vec ninja;
        const auto ref_magnitude = norm(x);

        // iteration counter
        auto counter = 0u;
        while(true) {
            const auto residual = func.evaluate_residual(x);

            if(0 == counter) ninja = solve(func.evaluate_jacobian(x), residual);
            else {
                // clear temporary factor container
                alpha.clear();
                alpha.reserve(hist_ninja.size());
                // commit current residual
                hist_residual.back() += residual;
                // commit current factor after obtaining ninja and residual
                hist_factor.emplace_back(dot(hist_ninja.back(), hist_residual.back()));
                // copy current residual to ninja
                ninja = residual;
                // perform two-step recursive loop
                // right side loop
                for(auto J = static_cast<int>(hist_factor.size()) - 1; J >= 0; --J) {
                    // compute and commit alpha
                    alpha.emplace_back(dot(hist_ninja[J], ninja) / hist_factor[J]);
                    // update ninja
                    ninja -= alpha.back() * hist_residual[J];
                }
                // apply the Hessian from the factorization in the first iteration
                ninja *= dot(hist_residual.back(), hist_ninja.back()) / dot(hist_residual.back(), hist_residual.back());
                // left side loop
                for(size_t I = 0, J = hist_factor.size() - 1; I < hist_factor.size(); ++I, --J) ninja += (alpha[J] - dot(hist_residual[I], ninja) / hist_factor[I]) * hist_ninja[I];
            }

            // commit current displacement increment
            hist_ninja.emplace_back(-ninja);
            hist_residual.emplace_back(-residual);

            const auto error = norm(ninja);
            const auto ref_error = error / ref_magnitude;
            suanpan_debug("LBFGS local iteration error: %.5E.\n", ref_error);
            if(error <= abs_tol && ref_error <= rel_tol) return SUANPAN_SUCCESS;
            if(++counter > max_iteration) return SUANPAN_FAIL;

            // check if the maximum record number is hit (LBFGS)
            if(counter > max_hist) {
                hist_ninja.pop_front();
                hist_residual.pop_front();
                hist_factor.pop_front();
            }

            x -= ninja;
        }
    }
};

class Quadratic {
public:
    [[nodiscard]] vec evaluate_residual(const vec& x) const { return square(x) - 1.; }

    [[nodiscard]] mat evaluate_jacobian(const vec& x) const { return diagmat(2. * x); }
};

#endif

//! @}
