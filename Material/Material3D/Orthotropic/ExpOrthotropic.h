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
 * @class ExpOrthotropic
 * @brief The ExpOrthotropic class.
 *
 * @author tlc
 * @date 20/02/2019
 * @version 0.1.0
 * @file ExpOrthotropic.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef EXPORTHOTROPIC_H
#define EXPORTHOTROPIC_H

#include "NonlinearOrthotropic.h"

struct DataExpOrthotropic {
    const double a, b;
};

class ExpOrthotropic : protected DataExpOrthotropic, public NonlinearOrthotropic {
    [[nodiscard]] double compute_k(const double p_strain) const override { return 1. + a - a * exp(-b * p_strain); }

    [[nodiscard]] double compute_dk(const double p_strain) const override { return a * b * exp(-b * p_strain); }

public:
    ExpOrthotropic(const unsigned T, const OrthotropicType TP, vec&& E, vec&& V, vec&& S, const double A, const double B, const double R)
        : DataExpOrthotropic{A, B}
        , NonlinearOrthotropic(T, TP, std::move(E), std::move(V), std::move(S), R) {}
};

class ExpHoffman final : public ExpOrthotropic {
public:
    ExpHoffman(const unsigned T, vec&& E, vec&& V, vec&& S, const double A, const double B, const double R)
        : ExpOrthotropic(T, OrthotropicType::Hoffman, std::move(E), std::move(V), std::move(S), A, B, R) {}

    unique_ptr<Material> get_copy() override { return std::make_unique<ExpHoffman>(*this); }
};

class ExpTsaiWu final : public ExpOrthotropic {
public:
    ExpTsaiWu(const unsigned T, vec&& E, vec&& V, vec&& S, const double A, const double B, const double R)
        : ExpOrthotropic(T, OrthotropicType::TsaiWu, std::move(E), std::move(V), std::move(S), A, B, R) {}

    unique_ptr<Material> get_copy() override { return std::make_unique<ExpTsaiWu>(*this); }
};

#endif

//! @}
