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
 * @class IntegrationPlan
 * @brief An IntegrationPlan class.
 *
 * This is a simple class provides users the integration points coordinates and
 * corresponding weights.
 *
 * Invoke the object with following definition:
 *
 *     IntegrationPlan TEST(intDimension, intOrder, intType);
 *
 * `intDimension` is the dimension of the problem.
 * `intOrder` is the integration order.
 * `intType` is the integration type. 1 is Gauss and 2 is Lobatto.
 *
 * The data is stored in `double** int_pts`.
 * The data can be acquired by using method `get_integration_scheme()`.
 * Single element could be accessed by `(int, int)` operator.
 * The data is arranged as follows.
 *
 *     |-------------------------------------|
 *     | POINT |   COORDINATES   |   WEIGHT  |
 *     |-------|-----|-----|-----|-----|-----|  <---
 *     |   1   | X_1 | Y_1 | ... |  WEIGHT_1 |    |
 *     |   2   | X_2 | Y_2 | ... |  WEIGHT_2 |    |
 *     |   3   | X_3 | Y_3 | ... |  WEIGHT_3 |    |
 *     |   .   |  .  |  .  |  .  |     .     |    |---ROWS=intOrder^intDimension
 *     |   .   |  .  |  .  |  .  |     .     |    |
 *     |   .   |  .  |  .  |  .  |     .     |    |
 *     |   N   | X_N | Y_N | ... |  WEIGHT_N |    |
 *     |-------------------------------------|  <---
 *             ^                             ^
 *             |-----------------------------|
 *                            |
 *                   COLS=intDimension+1
 *
 * @author tlc
 * @date 05/09/2017
 * @version 0.1.0
 * @file IntegrationPlan.h
 * @addtogroup Utility
 * @{
 */

#ifndef INTEGRATIONPLAN_H
#define INTEGRATIONPLAN_H

#include <armadillo/armadillo>

enum class IntegrationType {
    GAUSS,
    HERMITE,
    CHEBYSHEV,
    LOBATTO,
    RADAU,
    LAGUERRE,
    IRONS,
    TRIANGLE
};

class IntegrationPlan final {
    arma::mat int_pts;

    template<IntegrationType> void generate(unsigned, unsigned) { throw std::invalid_argument("not supported"); }

public:
    const unsigned n_rows = 0;
    const unsigned n_cols = 0;

    explicit IntegrationPlan(unsigned, unsigned, IntegrationType);

    [[nodiscard]] const arma::mat& get_data() const;

    double operator()(unsigned, unsigned) const;

    void print() const;
};

template<> void IntegrationPlan::generate<IntegrationType::GAUSS>(unsigned, unsigned);
template<> void IntegrationPlan::generate<IntegrationType::HERMITE>(unsigned, unsigned);
template<> void IntegrationPlan::generate<IntegrationType::CHEBYSHEV>(unsigned, unsigned);
template<> void IntegrationPlan::generate<IntegrationType::LOBATTO>(unsigned, unsigned);
template<> void IntegrationPlan::generate<IntegrationType::RADAU>(unsigned, unsigned);
template<> void IntegrationPlan::generate<IntegrationType::LAGUERRE>(unsigned, unsigned);
template<> void IntegrationPlan::generate<IntegrationType::IRONS>(unsigned, unsigned);
template<> void IntegrationPlan::generate<IntegrationType::TRIANGLE>(unsigned, unsigned);

#endif

//! @}
