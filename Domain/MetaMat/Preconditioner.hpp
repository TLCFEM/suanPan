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
 * @class Preconditioner
 * @brief A Preconditioner class.
 *
 * @author tlc
 * @date 21/07/2022
 * @version 0.1.0
 * @file Preconditioner.hpp
 * @addtogroup Preconditioner
 * @{
 */

#ifndef PRECONDITIONER_HPP
#define PRECONDITIONER_HPP

#include <suanPan.h>

template<sp_d data_t> class Preconditioner {
public:
    Preconditioner() = default;
    Preconditioner(const Preconditioner&) = default;
    Preconditioner(Preconditioner&&) noexcept = default;
    Preconditioner& operator=(const Preconditioner&) = default;
    Preconditioner& operator=(Preconditioner&&) noexcept = default;
    virtual ~Preconditioner() = default;

    virtual int init() { return SUANPAN_SUCCESS; }

    [[nodiscard]] virtual Col<data_t> apply(const Col<data_t>&) = 0;
};

template<sp_d data_t> class UnityPreconditioner final : public Preconditioner<data_t> {
public:
    [[nodiscard]] Col<data_t> apply(const Col<data_t>&) override;
};

template<sp_d data_t> Col<data_t> UnityPreconditioner<data_t>::apply(const Col<data_t>& in) { return in; }

#endif

//! @}
