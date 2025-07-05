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
 * @class CustomOrthotropic
 * @brief The CustomOrthotropic class.
 *
 * @author tlc
 * @date 16/01/2023
 * @version 0.2.0
 * @file CustomOrthotropic.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef CUSTOMORTHOTROPIC_H
#define CUSTOMORTHOTROPIC_H

#include "NonlinearOrthotropic.h"

#include <Toolbox/Expression.h>
#include <Toolbox/ResourceHolder.h>

class CustomOrthotropic : public NonlinearOrthotropic {
    [[nodiscard]] double compute_k(double) const override;
    [[nodiscard]] double compute_dk(double) const override;

    const unsigned k_tag;

    ResourceHolder<Expression> k_expression;

public:
    CustomOrthotropic(
        unsigned,        // tag
        OrthotropicType, // type
        vec&&,           // elastic modulus
        vec&&,           // poissons ratio
        vec&&,           // sigma
        unsigned,        // hardening function tag
        double = 0.      // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    void print() override;
};

class CustomHoffman final : public CustomOrthotropic {
public:
    CustomHoffman(const unsigned T, vec&& E, vec&& V, vec&& S, const unsigned K, const double R)
        : CustomOrthotropic(T, OrthotropicType::Hoffman, std::move(E), std::move(V), std::move(S), K, R) {}

    unique_ptr<Material> get_copy() override { return std::make_unique<CustomHoffman>(*this); }
};

class CustomTsaiWu final : public CustomOrthotropic {
public:
    CustomTsaiWu(const unsigned T, vec&& E, vec&& V, vec&& S, const unsigned K, const double R)
        : CustomOrthotropic(T, OrthotropicType::TsaiWu, std::move(E), std::move(V), std::move(S), K, R) {}

    unique_ptr<Material> get_copy() override { return std::make_unique<CustomTsaiWu>(*this); }
};

#endif

//! @}
