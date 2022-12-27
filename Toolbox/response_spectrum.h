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
 * @fn response_spectrum
 * @brief A function to compute response spectrum.
 *
 * Reference:
 *     1. https://doi.org/10.1016/S0267-7261(05)80015-6
 *
 * @author tlc
 * @date 20/10/2022
 * @version 0.1.0
 * @file response_spectrum.h
 * @addtogroup Utility
 * @{
 */

#ifndef RESPONSE_SPECTRUM_H
#define RESPONSE_SPECTRUM_H

#include <Toolbox/utility.h>

template<sp_d T> class Oscillator {
    const T omega; // natural angular frequency
    const T zeta;  // damping ratio

    const T alpha = omega * zeta;
    const T beta = omega * std::sqrt(1. - zeta * zeta);

    T gamma{0.};
    T a{0.}, b{0.}, c{0.};

    T amplitude(const Col<T>& data) { return std::max(std::abs(data.max()), std::abs(data.min())); }

    void compute_parameter(const T interval) {
        const auto exp_term = std::exp(-alpha * interval);

        a = exp_term * std::sin(beta * interval) / beta;
        b = 2. * exp_term * std::cos(beta * interval);
        c = exp_term * exp_term;

        gamma = (1. - b + c) / a / interval / omega / omega;
    }

    std::tuple<Col<T>, Col<T>, Col<T>> populate_response(const T interval, const Col<T>& motion) {
        this->compute_parameter(interval);

        Col<T> displacement(motion.n_elem, fill::none);
        displacement(0) = T(0);
        displacement(1) = b * displacement(0) - motion(0);

        for(auto I = 2llu, J = 1llu, K = 0llu; I < motion.n_elem; ++I, ++J, ++K) displacement(I) = b * displacement(J) - c * displacement(K) - motion(J);

        const auto n_elem = motion.n_elem - 1llu;

        Col<T> velocity(motion.n_elem, fill::none);
        velocity(0) = T(0);
        velocity.tail(n_elem) = arma::diff(displacement);

        Col<T> acceleration(motion.n_elem, fill::none);
        acceleration(0) = T(0);
        acceleration.tail(n_elem) = arma::diff(velocity);

        return std::make_tuple(displacement, velocity, acceleration);
    }

public:
    Oscillator(const T O, const T Z)
        : omega(O)
        , zeta(Z) {}

    Mat<T> compute_response(const T interval, const Col<T>& motion) {
        Col<T> displacement, velocity, acceleration;
        std::tie(displacement, velocity, acceleration) = this->populate_response(interval, motion);

        const auto factor = gamma * a;

        displacement *= factor * interval;
        velocity *= factor;
        acceleration *= factor / interval;

        return arma::join_rows(displacement, velocity, acceleration);
    }

    Col<T> compute_maximum_response(const T interval, const Col<T>& motion) {
        Col<T> displacement, velocity, acceleration;
        std::tie(displacement, velocity, acceleration) = this->populate_response(interval, motion);

        const auto factor = gamma * a;

        const auto max_u = this->amplitude(displacement) * factor * interval;
        const auto max_v = this->amplitude(velocity) * factor;
        const auto max_a = this->amplitude(acceleration * factor / interval + motion);

        return {max_u, max_v, max_a};
    }
};

/**
 * \brief compute response spectrum of the given ground motion
 * \param damping_ratio damping ratio
 * \param interval sampling interval of the target ground motion
 * \param motion target ground motion stored in one column
 * \param period periods where response spectrum needs to be computed
 * \return response spectrum stored in four columns (period, displacement, velocity, acceleration)
 */
template<sp_d T> Mat<T> response_spectrum(const T damping_ratio, const T interval, const Col<T>& motion, const Col<T>& period) {
    Mat<T> spectrum(3, period.n_elem, fill::none);

    suanpan_for(0llu, period.n_elem, [&](const uword I) {
        if(!suanpan::approx_equal(period(I), 0., 10000)) [[likely]] spectrum.col(I) = Oscillator(datum::tau / period(I), damping_ratio).compute_maximum_response(interval, motion);
        else [[unlikely]]
        {
            spectrum.col(I)(0) = spectrum.col(I)(1) = 0.;
            spectrum.col(I)(2) = std::max(std::abs(motion.max()), std::abs(motion.min()));
        }
    });

    arma::inplace_trans(spectrum);

    return arma::join_rows(period, spectrum);
}

/**
 * \brief compute response of a linear SDOF system
 * \param damping_ratio damping ratio of the system
 * \param interval time step size
 * \param freq frequency of the system (Hz)
 * \param motion acceleration record
 * \return response history
 */
template<sp_d T> Mat<T> sdof_response(const T damping_ratio, const T interval, const T freq, const Col<T>& motion) {
    Oscillator system(freq * datum::tau, damping_ratio);

    return arma::join_rows(interval * regspace(0., static_cast<double>(motion.n_elem) - 1.), system.compute_response(interval, motion));
}

#endif
