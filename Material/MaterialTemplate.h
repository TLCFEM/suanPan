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
 * @class MaterialTemplate
 * @brief A MaterialTemplate class.
 *
 * The Material class defines material models for FEM analysis.
 *
 * @author tlc
 * @date 05/09/2017
 * @version 0.1.1
 * @file MaterialTemplate.h
 * @addtogroup Material
 * @{
 */

#ifndef MATERIALTEMPLATE_H
#define MATERIALTEMPLATE_H

#include <Material/Material.h>

class MaterialTemplate final : public Material {
public:
    explicit MaterialTemplate(unsigned = 0);

    int initialize(const shared_ptr<DomainBase>&) override;

    [[nodiscard]] double get_parameter(ParameterType) const override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
