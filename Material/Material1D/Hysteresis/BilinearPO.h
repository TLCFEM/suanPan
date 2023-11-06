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
/**
 * @class BilinearPO
 * @brief A BilinearPO material class.
 * @author tlc
 * @date 14/04/2020
 * @version 0.1.0
 * @file BilinearPO.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef BILINEARPO_H
#define BILINEARPO_H

#include "PeakOriented.h"

struct DataBilinearPO {
    const double elastic_modulus;

    const double t_strain;
    const double t_hardening;
    const double c_strain;
    const double c_hardening;
};

class BilinearPO final : protected DataBilinearPO, public PeakOriented {
    const double t_stress = elastic_modulus * t_strain;
    const double c_stress = elastic_modulus * c_strain;

    [[nodiscard]] podarray<double> compute_compression_initial_reverse() const override;
    [[nodiscard]] podarray<double> compute_tension_initial_reverse() const override;
    [[nodiscard]] podarray<double> compute_tension_backbone(double) const override;
    [[nodiscard]] podarray<double> compute_compression_backbone(double) const override;

public:
    BilinearPO(
        int,    // tag
        double, // elastic modulus
        double, // tension yield strain
        double, // tension hardening ratio
        double, // compression yield strain
        double, // compression hardening ratio
        double  // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    void print() override;
};

#endif

//! @}
