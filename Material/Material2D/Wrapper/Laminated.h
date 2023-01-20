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
 * @class Laminated
 * @brief A Laminated material class.
 * @author tlc
 * @date 03/10/2017
 * @version 0.1.1
 * @file Laminated.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef LAMINATED_H
#define LAMINATED_H

#include <Material/Material2D/Material2D.h>

class Laminated final : public Material2D {
    const uvec mat_tag;

    std::vector<unique_ptr<Material>> mat_pool;

public:
    Laminated(unsigned, // tag
              uvec&&    // mat tag
    );
    Laminated(const Laminated&);
    Laminated(Laminated&&) = delete;
    Laminated& operator=(const Laminated&) = delete;
    Laminated& operator=(Laminated&&) = delete;
    ~Laminated() override = default;

    int initialize(const shared_ptr<DomainBase>&) override;

    unique_ptr<Material> get_copy() override;

    int update_trial_status(const vec&) override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;

    vector<vec> record(OutputType) override;

    void print() override;
};

#endif

//! @}
