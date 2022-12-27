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

// ReSharper disable CppClangTidyClangDiagnosticDocumentationUnknownCommand
#include "NLE3D01.h"
#include <Toolbox/tensorToolbox.h>

/**
 * \brief compute the derivatives of potential function
 * \param m_strain \tr(\varepsilon)
 * \param d_strain \dfrac{2}{3}\varepsilon_d:\varepsilon_d
 * \return a vector consists of partial derivatives
 */
vec NLE3D01::compute_derivative(const double m_strain, const double d_strain) {
    const auto eqv_strain_squared = std::max(datum::eps, d_strain);

    vec output(6);

    output(0) = bulk * m_strain;
    output(1) = factor_a * pow(eqv_strain_squared, .5 * m - .5);
    output(2) = bulk;
    output(3) = factor_b * pow(eqv_strain_squared, .5 * m - 1.5);
    output(4) = 0.;
    output(5) = 0.;

    return output;
}

NLE3D01::NLE3D01(const unsigned T, const double K, const double RE, const double RS, const double M, const double R)
    : DataNLE3D01{9 * K, fabs(RE), fabs(RS), std::max(0., std::min(1., M))}
    , IsotropicNonlinearElastic3D(T, R) {}

unique_ptr<Material> NLE3D01::get_copy() { return make_unique<NLE3D01>(*this); }
