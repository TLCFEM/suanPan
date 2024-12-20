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
 * @class ComplexHysteresis
 * @brief A ComplexHysteresis material class.
 *
 * @author tlc
 * @date 09/11/2023
 * @version 0.1.0
 * @file ComplexHysteresis.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef COMPLEXHYSTERESIS_H
#define COMPLEXHYSTERESIS_H

#include <Material/Material1D/Material1D.h>

class ComplexHysteresis : public Material1D {
    [[nodiscard]] virtual podarray<double> compute_compression_backbone(double) = 0;
    [[nodiscard]] virtual podarray<double> compute_tension_backbone(double) = 0;
    [[nodiscard]] virtual podarray<double> compute_compression_unload(double) = 0;
    [[nodiscard]] virtual podarray<double> compute_tension_unload(double) = 0;
    [[nodiscard]] virtual podarray<double> compute_compression_reload(double) = 0;
    [[nodiscard]] virtual podarray<double> compute_tension_reload(double) = 0;
    [[nodiscard]] virtual podarray<double> compute_compression_subunload(double) = 0;
    [[nodiscard]] virtual podarray<double> compute_tension_subunload(double) = 0;

    virtual void update_compression_unload(double) = 0;
    virtual void update_tension_unload(double) = 0;

protected:
    enum class Status {
        NONE,
        CBACKBONE,
        TBACKBONE,
        CUNLOAD,
        TUNLOAD,
        CSUBUNLOAD,
        TSUBUNLOAD,
        CRELOAD,
        TRELOAD,
        CTRANS,
        TTRANS
    };

    Status trial_load_status = Status::NONE, current_load_status = Status::NONE;

    const double elastic_modulus;

public:
    ComplexHysteresis(
        unsigned,   // tag
        double,     // elastic modulus
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
