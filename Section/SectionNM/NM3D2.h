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
 * @class NM3D2
 * @brief A NM3D2 class.
 * @author tlc
 * @date 26/11/2021
 * @version 0.1.0
 * @file NM3D2.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef NM3D2_H
#define NM3D2_H

#include <Section/SectionNM/NonlinearNM3D.h>

class NM3D2 final : public NonlinearNM3D {
    const mat para_set;
    const vec yield_force;
    const double c, h, k;

    [[nodiscard]] double evaluate(double, double, double, const mat&) const;
    [[nodiscard]] static vec differentiate(const mat&, uword, uword);

    [[nodiscard]] double compute_h(double) const override;
    [[nodiscard]] double compute_dh(double) const override;

    [[nodiscard]] double compute_f(const vec&) const override;
    [[nodiscard]] vec compute_df(const vec&) const override;
    [[nodiscard]] mat compute_ddf(const vec&) const override;

public:
    NM3D2(unsigned, // tag
          double,   // EA
          double,   // EIS
          double,   // EIW
          double,   // NP
          double,   // MPS
          double,   // MPW
          double,   // c
          double,   // h
          double,   // k
          double,   // linear density
          mat&& = {});

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Section> get_copy() override;

    void print() override;
};

#endif

//! @}
