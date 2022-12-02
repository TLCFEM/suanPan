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
 * @class RelResidual
 * @brief The Converger class handles converger test to indicate if the
 * iteration
 * converges according to various rules.
 *
 * The class stores a pointer `factory` pointed to the Workroom and get
 * information from
 * this Workroom. The `tolerance` and `error` are stored independently so
 * that the
 * Workroom will not be modified.
 *
 * The class further provides a `print_flag` to indicate if the test
 * information should be
 * printed out.
 *
 * @author tlc
 * @date 27/08/2017
 * @version 0.1.0
 * @file RelResidual.h
 * @addtogroup Converger
 * @{
 */

#ifndef RELRESIDUAL_H
#define RELRESIDUAL_H

#include "Converger.h"

class DomainBase;

class RelResidual final : public Converger {
    double ref_residual = 0.;

public:
    explicit RelResidual(unsigned = 0, double = 1E-8, unsigned = 7, bool = false);

    unique_ptr<Converger> get_copy() override;

    bool is_converged(unsigned) override;
};

#endif

//! @}
