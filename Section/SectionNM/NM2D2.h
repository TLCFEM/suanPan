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
 * @class NM2D2
 * @brief A NM2D2 class.
 * @author tlc
 * @date 30/11/2021
 * @version 0.1.0
 * @file NM2D2.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef NM2D2_H
#define NM2D2_H

#include "LinearHardeningNM.h"
#include "SurfaceNM2D.h"

class NM2D2 final : protected SurfaceNM2D, public LinearHardeningNM {
protected:
    [[nodiscard]] double compute_f(const vec&, const vec&) const override;
    [[nodiscard]] vec compute_df(const vec&, const vec&) const override;
    [[nodiscard]] mat compute_ddf(const vec&, const vec&) const override;

public:
    NM2D2(
        unsigned, // tag
        double,   // EA
        double,   // EIS
        double,   // NP
        double,   // MP
        double,   // c
        double,   // h
        double,   // k
        double,   // linear density
        mat&& = {}
    );

    unique_ptr<Section> get_copy() override;
};

#endif

//! @}
