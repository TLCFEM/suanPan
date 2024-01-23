/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @fn shape
 * @brief Provides common shape functions.
 *
 * @author tlc
 * @date 05/01/2022
 * @version 0.4.0
 * @file shape.h
 * @addtogroup Utility
 * @{
 */

#ifndef SHAPE_H
#define SHAPE_H

#include <suanPan.h>

namespace interpolation {
    template<typename T> Row<T> linear(const T X, const T Y) { return {1., X, Y, X * Y}; }

    template<typename T> Row<T> linear(const T X, const T Y, const T Z) { return {1., X, Y, Z, X * Y, Y * Z, Z * X}; }

    template<typename T> Row<T> linear(const Col<T>& C) { return C.size() == 2 ? linear(C(0), C(1)) : linear(C(0), C(1), C(2)); }

    template<typename T> Row<T> quadratic(const T X, const T Y) { return {1., X, Y, X * X, X * Y, Y * Y, X * X * Y, X * Y * Y, X * X * Y * Y}; }

    template<typename T> Row<T> quadratic(const Col<T>& C) { return quadratic(C(0), C(1)); }

    template<typename T> Row<T> cubic(const T X, const T Y) { return {1., X, Y, X * X, X * Y, Y * Y, X * X * X, X * X * Y, X * Y * Y, Y * Y * Y}; }

    template<typename T> Row<T> cubic(const Col<T>& C) { return cubic(C(0), C(1)); }
} // namespace interpolation

namespace area {
    template<typename T> T triangle(const Mat<T>& EC);

    template<typename T> T shoelace(const Mat<T>& C);
} // namespace area

namespace shape {
    template<typename T> Mat<T> truss(T int_pts, unsigned order = 0, unsigned num_node = 2);
    template<typename T> Col<T> beam(T int_pts, unsigned order, double length);
    template<typename T> Mat<T> triangle(const Col<T>& int_pts, unsigned order);
    template<typename T> Mat<T> quad(const Mat<T>& int_pts, unsigned order, unsigned num_node = 4);
    template<typename T> Mat<T> cube(const Mat<T>& int_pts, unsigned order, unsigned num_node = 8);

    namespace plate {
        template<typename T> Mat<T> triangle(const Col<T>& int_pts, unsigned order, unsigned num_node, const Mat<T>& nodes);
        template<typename T> Mat<T> quad(const Col<T>& int_pts, unsigned order, unsigned num_node = 4);
    } // namespace plate

    template<typename T> Mat<T> stress(T X, T Y, unsigned S);

    template<typename T> Mat<T> stress(const Col<T>& C, unsigned S);
    template<typename T> Mat<T> stress5(const Col<T>& C);
    template<typename T> Mat<T> stress7(const Col<T>& C);
    template<typename T> Mat<T> stress9(const Col<T>& C);
    template<typename T> Mat<T> stress11(const Col<T>& C);

    template<typename T> Mat<T> stress5(T X, T Y);
    template<typename T> Mat<T> stress7(T X, T Y);
    template<typename T> Mat<T> stress9(T X, T Y);
    template<typename T> Mat<T> stress11(T X, T Y);

    inline mat stress5(const vec& C);
    inline mat stress7(const vec& C);
    inline mat stress9(const vec& C);
    inline mat stress11(const vec& C);

    template<typename T> Mat<T> strain(T X, T Y, T V, unsigned S);

    template<typename T> Mat<T> strain(const Col<T>& C, T V, unsigned S);
    template<typename T> Mat<T> strain5(T X, T Y, T V);
    template<typename T> Mat<T> strain7(T X, T Y, T V);
    template<typename T> Mat<T> strain9(T X, T Y, T V);
    template<typename T> Mat<T> strain11(T X, T Y, T V);

    template<typename T> Mat<T> strain5(const Col<T>& C, T V);
    template<typename T> Mat<T> strain7(const Col<T>& C, T V);
    template<typename T> Mat<T> strain9(const Col<T>& C, T V);
    template<typename T> Mat<T> strain11(const Col<T>& C, T V);

    inline mat strain5(const vec& C, double V);
    inline mat strain7(const vec& C, double V);
    inline mat strain9(const vec& C, double V);
    inline mat strain11(const vec& C, double V);

    template<typename T> Mat<T> linear_stress(T X, T Y);
} // namespace shape

template<typename T> T area::triangle(const Mat<T>& EC) { return .5 * (EC(0, 0) * (EC(1, 1) - EC(2, 1)) + EC(1, 0) * (EC(2, 1) - EC(0, 1)) + EC(2, 0) * (EC(0, 1) - EC(1, 1))); }

template<typename T> T area::shoelace(const Mat<T>& C) {
    suanpan_assert([&] { if(2 != C.n_cols) throw invalid_argument("need two columns"); });

    const auto S = C.n_rows;
    Mat<T> E = arma::resize(C, S + 1, C.n_cols);

    E.tail_rows(1) = C.head_rows(1);

    const vec& X = E.col(0);
    const vec& Y = E.col(1);

    return .5 * fabs(arma::dot(X.head(S), Y.tail(S)) - arma::dot(X.tail(S), Y.head(S)));
}

template<typename T> Mat<T> shape::truss(const T int_pts, const unsigned order, const unsigned num_node) {
    Mat<T> N(1, num_node);

    if(const auto& X = int_pts; num_node == 2) {
        if(order == 0) {
            N(0, 0) = 1. - X;
            N(0, 1) = 1. + X;
        }
        else N(0, 0) = -(N(0, 1) = 1.);
        N *= .5;
    }
    else if(num_node == 3) {
        if(order == 0) {
            const auto XX = X * X;
            N(0, 0) = .5 * (XX - X);
            N(0, 1) = 1. - XX;
            N(0, 2) = .5 * (XX + X);
        }
        else {
            N(0, 0) = X - .5;
            N(0, 1) = -2. * X;
            N(0, 2) = X + .5;
        }
    }

    return N;
}

template<typename T> Col<T> shape::beam(const T int_pts, const unsigned order, const double length) {
    Col<T> N(4);

    const auto XP = 1. + int_pts;
    const auto XM = 1. - int_pts;
    const auto XPP = XP * XP;
    const auto XMM = XM * XM;

    if(order == 0) {
        N(0) = 2. * XMM * (XP + 1.);
        N(1) = length * XMM * XP;
        N(2) = 2. * XPP * (XM + 1.);
        N(3) = length * XM * XPP;
    }
    else if(order == 1) {
        N(0) = -6. * XP * XM;
        N(1) = length * XM * (3. * int_pts + 1.);
        N(2) = 6. * XP * XM;
        N(3) = length * XP * (3. * int_pts - 1.);
    }

    N *= .125;

    return N;
}

/**
 * \brief compute the shape function or its derivative of six node triangle in global coordinate system
 * \tparam T double/float
 * \param int_pts global coordinates of integration points
 * \param order shape function or its derivative
 * \return shape function or its derivative
 */
template<typename T> Mat<T> shape::triangle(const Col<T>& int_pts, const unsigned order) {
    Mat<T> N;

    if(order != 0 && order != 1) throw invalid_argument("order needs to be either 0 or 1");

    N.zeros(order + 1llu, 6);

    if(const auto &X = int_pts(0), &Y = int_pts(1); order == 0) {
        N(0, 0) = 1.;
        N(0, 1) = X;
        N(0, 2) = Y;
        N(0, 3) = X * Y;
        N(0, 4) = X * X;
        N(0, 5) = Y * Y;
    }
    else if(order == 1) {
        N(0, 1) = N(1, 2) = 1.;
        N(0, 4) = 2. * (N(1, 3) = X);
        N(1, 5) = 2. * (N(0, 3) = Y);
    }

    return N;
}

template<typename T> Mat<T> shape::quad(const Mat<T>& int_pts, const unsigned order, const unsigned num_node) {
    Mat<T> N;

    if(order != 0 && order != 1) throw invalid_argument("order needs to be either 0 or 1");
    if(num_node < 4 || num_node > 8) throw invalid_argument("number of nodes must between 4 and 8");

    N.zeros(order + 1llu, num_node);

    const auto& X = int_pts(0);
    const auto& Y = int_pts(1);

    if(const auto XP = 1. + X, XM = 1. - X, YP = 1. + Y, YM = 1. - Y; 8 == num_node) {
        if(const auto XX = X * X, YY = Y * Y, XY = X * Y; 0 == order) {
            N(0, 7) = .5 * XM * (1. - YY);
            N(0, 6) = .5 * (1. - XX) * YP;
            N(0, 5) = .5 * XP * (1. - YY);
            N(0, 4) = .5 * (1. - XX) * YM;
            N(0, 0) = .25 * XM * YM - .5 * (N(0, 4) + N(0, 7));
            N(0, 1) = .25 * XP * YM - .5 * (N(0, 4) + N(0, 5));
            N(0, 2) = .25 * XP * YP - .5 * (N(0, 5) + N(0, 6));
            N(0, 3) = .25 * XM * YP - .5 * (N(0, 6) + N(0, 7));
        }
        else if(1 == order) {
            const auto X2 = .5 * X;
            const auto Y2 = .5 * Y;
            const auto X4 = .25 * X;
            const auto Y4 = .25 * Y;
            const auto X24 = .25 * XX;
            const auto Y24 = .25 * YY;
            const auto XY2 = .5 * XY;
            N(1, 7) = XY - Y;
            N(1, 6) = .5 - .5 * XX;
            N(1, 5) = -Y - XY;
            N(1, 4) = .5 * XX - .5;
            N(1, 3) = Y2 - X4 - XY2 + X24;
            N(1, 2) = X4 + Y2 + XY2 + X24;
            N(1, 1) = Y2 - X4 + XY2 - X24;
            N(1, 0) = X4 + Y2 - XY2 - X24;
            N(0, 7) = .5 * YY - .5;
            N(0, 6) = -X - XY;
            N(0, 5) = .5 - .5 * YY;
            N(0, 4) = XY - X;
            N(0, 3) = X2 - Y4 + XY2 - Y24;
            N(0, 2) = X2 + Y4 + XY2 + Y24;
            N(0, 1) = X2 - Y4 - XY2 + Y24;
            N(0, 0) = X2 + Y4 - XY2 - Y24;
        }
    }
    else {
        if(0 == order) {
            N(0, 3) = XM * YP;
            N(0, 2) = XP * YP;
            N(0, 1) = XP * YM;
            N(0, 0) = XM * YM;
        }
        else if(1 == order) {
            N(1, 1) = -(N(1, 2) = XP);
            N(1, 0) = -(N(1, 3) = XM);
            N(0, 3) = -(N(0, 2) = YP);
            N(0, 0) = -(N(0, 1) = YM);
        }
        N *= .25;
        if(5 <= num_node) {
            if(0 == order) {
                N(0, 4) = .25 * (1. - X * X) * (1. - Y);
                N(0, 0) -= .5 * N(0, 4);
                N(0, 1) -= .5 * N(0, 4);
            }
            else {
                N(0, 4) = -.5 * X * (1. - Y);
                N(1, 4) = -.25 * (1. - X * X);
                N(0, 0) -= .5 * N(0, 4);
                N(0, 1) -= .5 * N(0, 4);
                N(1, 0) -= .5 * N(1, 4);
                N(1, 1) -= .5 * N(1, 4);
            }
        }
        if(6 <= num_node) {
            if(0 == order) {
                N(0, 5) = .25 * (1. - Y * Y) * (1. + X);
                N(0, 1) -= .5 * N(0, 5);
                N(0, 2) -= .5 * N(0, 5);
            }
            else {
                N(0, 5) = .25 * (1. - Y * Y);
                N(1, 5) = -.5 * Y * (1. + X);
                N(0, 1) -= .5 * N(0, 5);
                N(0, 2) -= .5 * N(0, 5);
                N(1, 1) -= .5 * N(1, 5);
                N(1, 2) -= .5 * N(1, 5);
            }
        }
        if(7 <= num_node) {
            if(0 == order) {
                N(0, 6) = .25 * (1. - X * X) * (1. + Y);
                N(0, 2) -= .5 * N(0, 6);
                N(0, 3) -= .5 * N(0, 6);
            }
            else {
                N(0, 6) = -.5 * X * (1. + Y);
                N(1, 6) = .25 * (1. - X * X);
                N(0, 2) -= .5 * N(0, 6);
                N(0, 3) -= .5 * N(0, 6);
                N(1, 2) -= .5 * N(1, 6);
                N(1, 3) -= .5 * N(1, 6);
            }
        }
    }

    return N;
}

template<typename T> Mat<T> shape::cube(const Mat<T>& int_pts, const unsigned order, const unsigned num_node) {
    Mat<T> N;

    if(order == 0) N.zeros(1, num_node);
    else if(order == 1) N.zeros(3, num_node);
    else throw invalid_argument("order needs to be either 0 or 1");

    const auto& X = int_pts(0);
    const auto& Y = int_pts(1);
    const auto& Z = int_pts(2);

    const auto XP = 1. + X;
    const auto XM = 1. - X;
    const auto YP = 1. + Y;
    const auto YM = 1. - Y;
    const auto ZP = 1. + Z;
    const auto ZM = 1. - Z;

    if(num_node == 8) {
        if(order == 0) {
            N(0, 0) = XM * YM * ZM;
            N(0, 1) = XP * YM * ZM;
            N(0, 2) = XP * YP * ZM;
            N(0, 3) = XM * YP * ZM;
            N(0, 4) = XM * YM * ZP;
            N(0, 5) = XP * YM * ZP;
            N(0, 6) = XP * YP * ZP;
            N(0, 7) = XM * YP * ZP;
        }
        else if(order == 1) {
            N(0, 0) = -YM * ZM;
            N(0, 1) = YM * ZM;
            N(0, 2) = YP * ZM;
            N(0, 3) = -YP * ZM;
            N(0, 4) = -YM * ZP;
            N(0, 5) = YM * ZP;
            N(0, 6) = YP * ZP;
            N(0, 7) = -YP * ZP;
            N(1, 0) = -XM * ZM;
            N(1, 1) = -XP * ZM;
            N(1, 2) = XP * ZM;
            N(1, 3) = XM * ZM;
            N(1, 4) = -XM * ZP;
            N(1, 5) = -XP * ZP;
            N(1, 6) = XP * ZP;
            N(1, 7) = XM * ZP;
            N(2, 0) = -XM * YM;
            N(2, 1) = -XP * YM;
            N(2, 2) = -XP * YP;
            N(2, 3) = -XM * YP;
            N(2, 4) = XM * YM;
            N(2, 5) = XP * YM;
            N(2, 6) = XP * YP;
            N(2, 7) = XM * YP;
        }
        N *= .125;

        return N;
    }

    if(num_node == 20) {
        if(const auto XX = XP * XM, YY = YP * YM, ZZ = ZP * ZM; order == 0) {
            N(0, 0) = .125 * XM * YM * ZM * (-2. - X - Y - Z);
            N(0, 1) = .125 * XP * YM * ZM * (-2. + X - Y - Z);
            N(0, 2) = .125 * XP * YP * ZM * (-2. + X + Y - Z);
            N(0, 3) = .125 * XM * YP * ZM * (-2. - X + Y - Z);
            N(0, 4) = .125 * XM * YM * ZP * (-2. - X - Y + Z);
            N(0, 5) = .125 * XP * YM * ZP * (-2. + X - Y + Z);
            N(0, 6) = .125 * XP * YP * ZP * (-2. + X + Y + Z);
            N(0, 7) = .125 * XM * YP * ZP * (-2. - X + Y + Z);
            N(0, 8) = .25 * XX * YM * ZM;
            N(0, 9) = .25 * YY * XP * ZM;
            N(0, 10) = .25 * XX * YP * ZM;
            N(0, 11) = .25 * YY * XM * ZM;
            N(0, 12) = .25 * XX * YM * ZP;
            N(0, 13) = .25 * YY * XP * ZP;
            N(0, 14) = .25 * XX * YP * ZP;
            N(0, 15) = .25 * YY * XM * ZP;
            N(0, 16) = .25 * ZZ * XM * YM;
            N(0, 17) = .25 * ZZ * XP * YM;
            N(0, 18) = .25 * ZZ * XP * YP;
            N(0, 19) = .25 * ZZ * XM * YP;
        }
        else if(order == 1) {
            N(0, 0) = YM * ZM * (2. * X + Y + Z + 1.) * .125;
            N(0, 1) = YM * ZM * (2. * X - Y - Z - 1.) * .125;
            N(0, 2) = YP * ZM * (2. * X + Y - Z - 1.) * .125;
            N(0, 3) = YP * ZM * (2. * X - Y + Z + 1.) * .125;
            N(0, 4) = YM * ZP * (2. * X + Y - Z + 1.) * .125;
            N(0, 5) = YM * ZP * (2. * X - Y + Z - 1.) * .125;
            N(0, 6) = YP * ZP * (2. * X + Y + Z - 1.) * .125;
            N(0, 7) = YP * ZP * (2. * X - Y - Z + 1.) * .125;
            N(1, 0) = XM * ZM * (2. * Y + X + Z + 1.) * .125;
            N(1, 1) = XP * ZM * (2. * Y - X + Z + 1.) * .125;
            N(1, 2) = XP * ZM * (2. * Y + X - Z - 1.) * .125;
            N(1, 3) = XM * ZM * (2. * Y - X - Z - 1.) * .125;
            N(1, 4) = XM * ZP * (2. * Y + X - Z + 1.) * .125;
            N(1, 5) = XP * ZP * (2. * Y - X - Z + 1.) * .125;
            N(1, 6) = XP * ZP * (2. * Y + X + Z - 1.) * .125;
            N(1, 7) = XM * ZP * (2. * Y - X + Z - 1.) * .125;
            N(2, 0) = XM * YM * (2. * Z + X + Y + 1.) * .125;
            N(2, 1) = XP * YM * (2. * Z + Y - X + 1.) * .125;
            N(2, 2) = XP * YP * (2. * Z - X - Y + 1.) * .125;
            N(2, 3) = XM * YP * (2. * Z + X - Y + 1.) * .125;
            N(2, 4) = XM * YM * (2. * Z - X - Y - 1.) * .125;
            N(2, 5) = XP * YM * (2. * Z + X - Y - 1.) * .125;
            N(2, 6) = XP * YP * (2. * Z + X + Y - 1.) * .125;
            N(2, 7) = XM * YP * (2. * Z - X + Y - 1.) * .125;

            N(0, 8) = -X * YM * ZM * .5;
            N(0, 9) = YY * ZM * .25;
            N(0, 10) = -X * YP * ZM * .5;
            N(0, 11) = -YY * ZM * .25;
            N(1, 8) = -XX * ZM * .25;
            N(1, 9) = -Y * XP * ZM * .5;
            N(1, 10) = XX * ZM * .25;
            N(1, 11) = -Y * XM * ZM * .5;
            N(2, 8) = -XX * YM * .25;
            N(2, 9) = -XP * YY * .25;
            N(2, 10) = -XX * YP * .25;
            N(2, 11) = -XM * YY * .25;

            N(0, 12) = -X * YM * ZP * .5;
            N(0, 13) = YY * ZP * .25;
            N(0, 14) = -X * YP * ZP * .5;
            N(0, 15) = -YY * ZP * .25;
            N(1, 12) = -XX * ZP * .25;
            N(1, 13) = -Y * XP * ZP * .5;
            N(1, 14) = XX * ZP * .25;
            N(1, 15) = -Y * XM * ZP * .5;
            N(2, 12) = XX * YM * .25;
            N(2, 13) = XP * YY * .25;
            N(2, 14) = XX * YP * .25;
            N(2, 15) = XM * YY * .25;

            N(0, 16) = -YM * ZZ * .25;
            N(0, 17) = YM * ZZ * .25;
            N(0, 18) = YP * ZZ * .25;
            N(0, 19) = -YP * ZZ * .25;
            N(1, 16) = -XM * ZZ * .25;
            N(1, 17) = -XP * ZZ * .25;
            N(1, 18) = XP * ZZ * .25;
            N(1, 19) = XM * ZZ * .25;
            N(2, 16) = -Z * XM * YM * .5;
            N(2, 17) = -Z * XP * YM * .5;
            N(2, 18) = -Z * XP * YP * .5;
            N(2, 19) = -Z * XM * YP * .5;
        }

        return N;
    }

    throw invalid_argument("not supported");
}

template<typename T> Mat<T> shape::stress(const T X, const T Y, const unsigned S) {
    Mat<T> N = zeros(3, S);

    for(auto I = 0; I < 3; ++I) N(I, I) = 1.;

    if(S >= 5) {
        N(0, 4) = Y;
        N(1, 3) = X;
        if(S >= 7) {
            N(2, 5) = -(N(0, 6) = X);
            N(2, 6) = -(N(1, 5) = Y);
            if(S >= 9) {
                const auto X2 = X * X;
                const auto Y2 = Y * Y;
                const auto XY = X * Y;
                N(1, 7) = N(0, 8) = 2. * XY;
                N(2, 7) = -X2;
                N(2, 8) = -Y2;
                if(S == 11) {
                    N(1, 9) = 2. * X2 + (N(1, 10) = -Y2);
                    N(0, 10) = 2. * Y2 + (N(0, 9) = -X2);
                    N(2, 10) = N(2, 9) = 2. * XY;
                }
            }
        }
    }

    return N;
}

template<typename T> Mat<T> shape::stress(const Col<T>& C, const unsigned S) { return stress(C(0), C(1), S); }

template<typename T> Mat<T> shape::stress5(const Col<T>& C) { return stress(C, 5); }

template<typename T> Mat<T> shape::stress7(const Col<T>& C) { return stress(C, 7); }

template<typename T> Mat<T> shape::stress9(const Col<T>& C) { return stress(C, 9); }

template<typename T> Mat<T> shape::stress11(const Col<T>& C) { return stress(C, 11); }

template<typename T> Mat<T> shape::stress5(const T X, const T Y) { return stress(X, Y, 5); }

template<typename T> Mat<T> shape::stress7(const T X, const T Y) { return stress(X, Y, 7); }

template<typename T> Mat<T> shape::stress9(const T X, const T Y) { return stress(X, Y, 9); }

template<typename T> Mat<T> shape::stress11(const T X, const T Y) { return stress(X, Y, 11); }

template<typename T> Mat<T> shape::strain(const T X, const T Y, const T V, const unsigned S) {
    Mat<T> N(3, S, fill::zeros);

    N(0, 0) = N(1, 1) = 1.;

    N(2, 2) = 2. + 2. * V;

    N(0, 1) = N(1, 0) = -V;

    if(S >= 5) {
        N(0, 3) = -V * (N(1, 3) = X);
        N(1, 4) = -V * (N(0, 4) = Y);
        if(S >= 7) {
            N(0, 5) = N(1, 4);
            N(0, 6) = N(1, 3);

            N(1, 5) = N(0, 4);
            N(1, 6) = N(0, 3);

            N(2, 5) = -X * N(2, 2);
            N(2, 6) = -Y * N(2, 2);
            if(S >= 9) {
                const auto X2 = X * X, Y2 = Y * Y, XY = X * Y;

                N(1, 8) = N(0, 7) = -V * (N(1, 7) = N(0, 8) = 2. * XY);

                N(2, 7) = -X2 * N(2, 2);
                N(2, 8) = -Y2 * N(2, 2);
                if(S == 11) {
                    N(0, 9) = V * Y2 - (2. * V + 1.) * X2;
                    N(1, 9) = (2. + V) * X2 - Y2;

                    N(0, 10) = (2. + V) * Y2 - X2;
                    N(1, 10) = V * X2 - (2. * V + 1.) * Y2;

                    N(2, 10) = N(2, 9) = 2. * XY * N(2, 2);
                }
            }
        }
    }

    return N;
}

template<typename T> Mat<T> shape::strain(const Col<T>& C, const T V, const unsigned S) { return strain(C(0), C(1), V, S); }

template<typename T> Mat<T> shape::strain5(const T X, const T Y, const T V) { return strain(X, Y, V, 5); }

template<typename T> Mat<T> shape::strain7(const T X, const T Y, const T V) { return strain(X, Y, V, 7); }

template<typename T> Mat<T> shape::strain9(const T X, const T Y, const T V) { return strain(X, Y, V, 9); }

template<typename T> Mat<T> shape::strain11(const T X, const T Y, const T V) { return strain(X, Y, V, 11); }

template<typename T> Mat<T> shape::strain5(const Col<T>& C, const T V) { return strain(C, V, 5); }

template<typename T> Mat<T> shape::strain7(const Col<T>& C, const T V) { return strain(C, V, 7); }

template<typename T> Mat<T> shape::strain9(const Col<T>& C, const T V) { return strain(C, V, 9); }

template<typename T> Mat<T> shape::strain11(const Col<T>& C, const T V) { return strain(C, V, 11); }

template<typename T> Mat<T> shape::linear_stress(T X, T Y) {
    Mat<T> N = eye(3, 9);
    N.cols(3, 5) = X * eye(3, 3);
    N.cols(6, 8) = Y * eye(3, 3);

    return N;
}

template<typename T> Mat<T> shape::plate::triangle(const Col<T>& int_pts, const unsigned order, const unsigned num_node, const Mat<T>& nodes) {
    suanpan_assert([&] { if(order > 1 || (num_node != 3 && num_node != 6) || (nodes.n_cols != 2) || nodes.n_rows < 3) throw invalid_argument("not supported"); });

    Mat<T> N;

    if(order == 0) N.zeros(1, 3llu * num_node);
    else if(order == 1) N.zeros(3, 3llu * num_node);
    else throw invalid_argument("order needs to be either 0 or 1");

    Mat<T> TEMP(3, 3);
    TEMP.row(0).fill(1.);
    TEMP.rows(1, 2) = nodes.rows(0, 2).t();

    const Mat<T> A = inv(TEMP);
    const Col<T> W = solve(TEMP, vec{1., int_pts(0), int_pts(1)});

    const vec L = std::initializer_list<double>{W(0), W(1), W(2), W(0), W(1)};
    const vec B = std::initializer_list<double>{A(0, 1), A(1, 1), A(2, 1), A(0, 1), A(1, 1)};
    const vec C = std::initializer_list<double>{A(0, 2), A(1, 2), A(2, 2), A(0, 2), A(1, 2)};

    const auto DA = A.cols(1, 2);

    if(3 == num_node) {
        if(0 == order) {
            auto IDX = 0;
            for(auto I = 0; I < 3; ++I) {
                const auto J = I + 1;
                const auto K = J + 1;
                N(IDX++) = L(I) * (1. - L(J) * L(J) - L(K) * L(K)) + L(I) * L(I) * (L(J) + L(K));
                N(IDX++) = B(J) * L(I) * L(K) * (L(I) + .5 * L(J)) - B(K) * L(I) * L(J) * (L(I) + .5 * L(K));
                N(IDX++) = C(J) * L(I) * L(K) * (L(I) + .5 * L(J)) - C(K) * L(I) * L(J) * (L(I) + .5 * L(K));
            }
        }
        else {
            mat DNDL(3, 3);
            auto IDX = 0;
            for(auto I = 0; I < 3; ++I) {
                const auto J = I + 1;
                const auto K = J + 1;

                DNDL(0, 0) = L(J) + L(K);
                DNDL(1, 1) = DNDL(2, 2) = -L(I);
                DNDL(0, 1) = DNDL(1, 0) = L(I) - L(J);
                DNDL(0, 2) = DNDL(2, 0) = L(I) - L(K);
                DNDL(1, 2) = DNDL(2, 1) = 0.;

                mat DD = DA.t() * DNDL * DA;

                N(0, IDX) = -2. * DD(0, 0);
                N(1, IDX) = -2. * DD(1, 1);
                N(2, IDX++) = -2. * (DD(0, 1) + DD(1, 0));

                DNDL(0, 0) = 4. * (B(J) * L(K) - B(K) * L(J));
                DNDL(1, 1) = DNDL(2, 2) = 0.;
                DNDL(0, 1) = DNDL(1, 0) = (B(J) - B(K)) * L(K) - 4 * B(K) * L(I);
                DNDL(0, 2) = DNDL(2, 0) = 4. * B(J) * L(I) + (B(J) - B(K)) * L(J);
                DNDL(1, 2) = DNDL(2, 1) = L(I) * (B(J) - B(K));

                DD = DA.t() * DNDL * DA;

                N(0, IDX) = -.5 * DD(0, 0);
                N(1, IDX) = -.5 * DD(1, 1);
                N(2, IDX++) = -.5 * (DD(0, 1) + DD(1, 0));

                DNDL(0, 0) = 4. * (C(J) * L(K) - C(K) * L(J));
                DNDL(0, 1) = DNDL(1, 0) = (C(J) - C(K)) * L(K) - 4 * C(K) * L(I);
                DNDL(0, 2) = DNDL(2, 0) = 4. * C(J) * L(I) + (C(J) - C(K)) * L(J);
                DNDL(1, 2) = DNDL(2, 1) = L(I) * (C(J) - C(K));

                DD = DA.t() * DNDL * DA;

                N(0, IDX) = -.5 * DD(0, 0);
                N(1, IDX) = -.5 * DD(1, 1);
                N(2, IDX++) = -.5 * (DD(0, 1) + DD(1, 0));
            }
        }
    }

    return N;
}

template<typename T> Mat<T> shape::plate::quad(const Col<T>& int_pts, const unsigned order, const unsigned num_node) {
    Mat<T> N;

    if(order == 0) N.zeros(1, 3llu * num_node);
    else if(order == 1) N.zeros(2, 6llu * num_node);
    else throw invalid_argument("order needs to be either 0 or 1");

    const auto& X = int_pts(0);
    const auto& Y = int_pts(1);

    return N;
}

mat shape::stress5(const vec& C) { return stress(C, 5); }

mat shape::stress7(const vec& C) { return stress(C, 7); }

mat shape::stress9(const vec& C) { return stress(C, 9); }

mat shape::stress11(const vec& C) { return stress(C, 11); }

mat shape::strain5(const vec& C, const double V) { return strain(C, V, 5); }

mat shape::strain7(const vec& C, const double V) { return strain(C, V, 7); }

mat shape::strain9(const vec& C, const double V) { return strain(C, V, 9); }

mat shape::strain11(const vec& C, const double V) { return strain(C, V, 11); }

#endif

//! @}
