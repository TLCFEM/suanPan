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
 * @class ExternalMaterial
 * @brief A ExternalMaterial class.
 *
 * The ExternalMaterial class enables communication across DLL boundary via C style interface.
 * The function has such a signature:
 * ```
 * void(*)(ExternalMaterialData* data, int* info);
 * ```
 *
 * On enter, `info` indicates what kind of operation shall the external model do.
 * The following values are used:
 * 0: acquire memory and initialisation
 * 1: deallocate memory
 * 2: update trial state with strain input only
 * 3: update trial state with strain and strain rate inputs
 * 4: commit trial state
 * 5: reset to current state
 * 6: clear to empty state
 * 7: validate constants
 *
 * On exit, `info` indicates if the operation succeeds.
 * The following value is used:
 * 0: success
 *
 * @author tlc
 * @date 15/01/2021
 * @version 0.1.0
 * @file ExternalMaterial.h
 * @addtogroup Material
 * @{
 */

#ifndef EXTERNALMATERIAL_H
#define EXTERNALMATERIAL_H

#include <Material/Material.h>
#include "ExternalMaterialData.h"

class ExternalMaterial final : public Material {
    static MaterialType get_type(const ExternalMaterialData&);

    using Interface = void (*)(ExternalMaterialData*, int*);

    std::vector<double> constant;

    Interface cooker;

    ExternalMaterialData data;

public:
    ExternalMaterial(unsigned,              // unique material tag
                     std::vector<double>&&, // parameter pool
                     void*                  // handler pointer
    );
    ExternalMaterial(const ExternalMaterial&);
    ExternalMaterial(ExternalMaterial&&) noexcept = delete;
    ExternalMaterial& operator=(const ExternalMaterial&) = delete;
    ExternalMaterial& operator=(ExternalMaterial&&) noexcept = delete;
    ~ExternalMaterial() override;

    bool validate();

    int initialize(const shared_ptr<DomainBase>&) override;

    void initialize_history(unsigned) override;
    void set_initial_history(const vec&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;
    int update_trial_status(const vec&, const vec&) override;

    int commit_status() override;
    int reset_status() override;
    int clear_status() override;

    std::vector<vec> record(OutputType) override;
};

#endif

//! @}
