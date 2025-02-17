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
 * @class MultilinearOO
 * @brief A MultilinearOO material class.
 * @author tlc
 * @date 23/08/2020
 * @version 0.1.0
 * @file MultilinearOO.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef MULTILINEAROO_H
#define MULTILINEAROO_H

#include <Material/Material1D/Hysteresis/OriginOriented.h>

struct DataMultilinearOO {
    const mat t_backbone;
    const mat c_backbone;
};

class MultilinearOO final : protected DataMultilinearOO, public OriginOriented {
    [[nodiscard]] pod2 compute_tension_backbone(double) const override;
    [[nodiscard]] pod2 compute_compression_backbone(double) const override;

public:
    MultilinearOO(int, mat&&, mat&&, double);

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
