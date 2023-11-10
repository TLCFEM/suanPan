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
 * @class Concrete21
 * @brief A Concrete21 material class.
 *
 * For plain concrete.
 * Only monotonic load is supported at the moment.
 * Cyclic behaviour does not converge in most cases.
 *
 * algorithm validated @ 05/02/2019 by tlc
 *
 * @author tlc
 * @date 04/02/2019
 * @version 1.0.0
 * @file Concrete21.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef CONCRETE21_H
#define CONCRETE21_H

#include <Material/Material1D/Concrete/ConcreteTsai.h>
#include <Material/Material2D/Material2D.h>

class Concrete21 final : public Material2D {
    ConcreteTsai concrete_major, concrete_minor;

    const double elastic_modulus;
    double shear_modulus = 0.;

public:
    Concrete21(
        unsigned,   // tag
        double,     // elastic modulus
        double,     // peak stress in negative
        double,     // crack stress in positive
        double,     // NC
        double,     // NT
        double,     // middle point
        double,     // peak strain in negative
        double,     // crack strain in positive
        double = 0. // density
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    void print() override;
};

#endif

//! @}
