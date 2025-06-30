/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
 * @class BilinearOrthotropic
 * @brief The BilinearOrthotropic class.
 *
 * @author tlc
 * @date 20/01/2019
 * @version 0.2.0
 * @file BilinearOrthotropic.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef BILINEARORTHOTROPIC_H
#define BILINEARORTHOTROPIC_H

#include "NonlinearOrthotropic.h"

struct DataBilinearOrthotropic {
    const double hardening_modulus;
};

class BilinearOrthotropic : protected DataBilinearOrthotropic, public NonlinearOrthotropic {
    [[nodiscard]] double compute_k(const double p_strain) const override { return hardening_modulus >= 0. || p_strain <= -1. / hardening_modulus ? 1. + p_strain * hardening_modulus : 0.; }

    [[nodiscard]] double compute_dk(const double p_strain) const override { return hardening_modulus >= 0. || p_strain <= -1. / hardening_modulus ? hardening_modulus : 0.; }

public:
    BilinearOrthotropic(const unsigned T, const OrthotropicType TP, vec&& E, vec&& V, vec&& S, const double H, const double R)
        : DataBilinearOrthotropic{H}
        , NonlinearOrthotropic(T, TP, std::move(E), std::move(V), std::move(S), R) {}
};

class BilinearHoffman : public BilinearOrthotropic {
public:
    BilinearHoffman(const unsigned T, vec&& E, vec&& V, vec&& S, const double H, const double R)
        : BilinearOrthotropic(T, OrthotropicType::Hoffman, std::move(E), std::move(V), std::move(S), H, R) {}

    unique_ptr<Material> get_copy() override { return std::make_unique<BilinearHoffman>(*this); }
};

class BilinearTsaiWu : public BilinearOrthotropic {
public:
    BilinearTsaiWu(const unsigned T, vec&& E, vec&& V, vec&& S, const double H, const double R)
        : BilinearOrthotropic(T, OrthotropicType::TsaiWu, std::move(E), std::move(V), std::move(S), H, R) {}

    unique_ptr<Material> get_copy() override { return std::make_unique<BilinearTsaiWu>(*this); }
};

#endif

//! @}
