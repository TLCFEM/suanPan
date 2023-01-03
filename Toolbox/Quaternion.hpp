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
 * @class Quaternion
 * @brief An Quaternion class.
 *
 * @author tlc
 * @date 05/09/2020
 * @version 0.1.0
 * @file Quaternion.hpp
 * @addtogroup Utility
 * @{
 */

#ifndef QUATERNION_H
#define QUATERNION_H

#include <Toolbox/tensorToolbox.h>

template<typename T> class Quaternion {
    T re;
    Col<T> im;

public:
    Quaternion();
    Quaternion(T, T, T, T);
    Quaternion(T, const Col<T>&);
    Quaternion(T, Col<T>&&);

    const T& real() const;
    const Col<T>& imag() const;

    T norm() const;
    Quaternion& normalise();

    [[nodiscard]] Quaternion inv() const;
    [[nodiscard]] Quaternion conj() const;

    Quaternion operator+(const Quaternion&) const;
    Quaternion& operator+=(const Quaternion&);
    Quaternion operator-(const Quaternion&) const;
    Quaternion& operator-=(const Quaternion&);
    Quaternion operator*(const Quaternion&) const;
    Quaternion& operator*=(const Quaternion&);
    Quaternion operator/(const Quaternion&) const;
    Quaternion& operator/=(const Quaternion&);

    void print() const;

    friend bool operator==(const Quaternion& A, const Quaternion& B) {
        constexpr int ulp = 10000;
        if(suanpan::approx_equal(A.re, B.re, ulp) && suanpan::approx_equal(A.im(0), B.im(0), ulp) && suanpan::approx_equal(A.im(1), B.im(1), ulp) && suanpan::approx_equal(A.im(2), B.im(2), ulp)) return true;

        if(suanpan::approx_equal(A.re, -B.re, ulp) && suanpan::approx_equal(A.im(0), -B.im(0), ulp) && suanpan::approx_equal(A.im(1), -B.im(1), ulp) && suanpan::approx_equal(A.im(2), -B.im(2), ulp)) return true;

        return false;
    }

    friend Quaternion operator-(const Quaternion& A) { return Quaternion(-A.re, -A.im); }

    friend Quaternion operator-(Quaternion&& A) {
        A.re = -A.re;
        A.im *= -1.;
        return std::forward<Quaternion>(A);
    }

    Mat<T> operator*(const Mat<T>&) const;

    Mat<T> to_mat() const;
    Col<T> to_pseudo() const;
};

template<typename T> Quaternion<T>::Quaternion()
    : re(T(0))
    , im(arma::zeros<Col<T>>(3)) {}

template<typename T> Quaternion<T>::Quaternion(const T R, const T I, const T J, const T K)
    : re(R)
    , im({I, J, K}) {}

template<typename T> Quaternion<T>::Quaternion(const T R, const Col<T>& I)
    : re(R)
    , im(I) {}

template<typename T> Quaternion<T>::Quaternion(const T R, Col<T>&& I)
    : re(R)
    , im(std::forward<Col<T>>(I)) {}

template<typename T> const T& Quaternion<T>::real() const { return re; }

template<typename T> const Col<T>& Quaternion<T>::imag() const { return im; }

template<typename T> T Quaternion<T>::norm() const { return re * re + arma::dot(im, im); }

template<typename T> Quaternion<T>& Quaternion<T>::normalise() {
    const auto magnitude = std::sqrt(norm());
    re /= magnitude;
    im /= magnitude;
    return *this;
}

template<typename T> Quaternion<T> Quaternion<T>::inv() const {
    const auto L = norm();

    return Quaternion<T>(re / L, -im / L);
}

template<typename T> Quaternion<T> Quaternion<T>::conj() const { return Quaternion<T>(re, -im); }

template<typename T> Quaternion<T> Quaternion<T>::operator+(const Quaternion& B) const {
    Quaternion<T> A = *this;

    return A += B;
}

template<typename T> Quaternion<T>& Quaternion<T>::operator+=(const Quaternion& B) {
    re += B.re;
    im += B.im;

    return *this;
}

template<typename T> Quaternion<T> Quaternion<T>::operator-(const Quaternion& B) const {
    Quaternion<T> A = *this;

    return A -= B;
}

template<typename T> Quaternion<T>& Quaternion<T>::operator-=(const Quaternion& B) {
    re -= B.re;
    im -= B.im;

    return *this;
}

template<typename T> Quaternion<T> Quaternion<T>::operator*(const Quaternion& B) const {
    Quaternion<T> A;

    A.re = re * B.re - arma::dot(im, B.im);
    A.im = re * B.im + B.re * im + arma::cross(im, B.im);

    return A;
}

template<typename T> Quaternion<T>& Quaternion<T>::operator*=(const Quaternion& B) { return *this = *this * B; }

template<typename T> Quaternion<T> Quaternion<T>::operator/(const Quaternion& B) const { return *this * B.inv(); }

template<typename T> Quaternion<T>& Quaternion<T>::operator/=(const Quaternion& B) { return *this = *this * B.inv(); }

template<typename T> void Quaternion<T>::print() const { sp_info("re: {:+0.6E} im: {:+0.6E} {:+0.6E} {:+0.6E}\n", re, im(0), im(1), im(2)); }

template<typename T> Mat<T> Quaternion<T>::operator*(const Mat<T>& I) const { return to_mat() * I; }

template<typename T> Mat<T> Quaternion<T>::to_mat() const { return 2. * re * transform::skew_symm(im) + 2. * im * im.t() + (re * re - arma::dot(im, im)) * eye(3, 3); }

template<typename T> Col<T> Quaternion<T>::to_pseudo() const {
    const auto norm_im = arma::norm(im);

    if(suanpan::approx_equal(norm_im, 0., 1000)) return Col<T>(3, fill::zeros);

    vec rotation = re < 0. ? -im : im;

    rotation *= 2. / norm_im * (norm_im < std::abs(re) ? std::asin(norm_im) : std::acos(std::abs(re)));

    return rotation;
}

#endif

//! @}
